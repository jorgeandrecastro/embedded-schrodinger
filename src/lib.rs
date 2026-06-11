// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2025 Andre
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

#![no_std]
#![forbid(unsafe_code)]
#![doc = include_str!("../README.md")]

use embedded_f32_sqrt::sqrt;

/// Solveur de l'équation de Schrödinger **indépendante du temps en 1D**.
///
/// La résolution est basée sur la **méthode de tir** (*shooting method*)
/// combinée à une dichotomie binaire de 30 itérations.
///
/// L'implémentation est entièrement `no_std`, sans `libm` ni liaison C.
/// Toute opération de racine carrée est déléguée à
/// [`embedded-f32-sqrt`](https://docs.rs/embedded-f32-sqrt), qui utilise
/// l'algorithme de Newton-Raphson sur `f32` IEEE 754.
///
/// # Paramètre générique
///
/// `N` : nombre de points de la grille spatiale (doit être ≥ 3).
///
/// # Exemple
///
/// ```rust
/// use embedded_schrodinger::SchrodingerSolver;
///
/// const N: usize = 100;
///
/// // Puits de potentiel infini (V = 0 à l'intérieur)
/// let potential = [0.0f32; N];
///
/// // Unités atomiques : ħ = 1, m = 1, dx = 0.01
/// let solver = SchrodingerSolver::new(potential, 0.01, 1.0, 1.0);
///
/// let mut psi = [0.0f32; N];
/// let energy = solver.find_eigenstate(0.1, 20.0, &mut psi)
///     .expect("Pas de solution trouvée dans l'intervalle donné");
///
/// // En unités atomiques, E₁ ≈ π²/2 ≈ 4.935 pour L = 1 (N·dx = 1.0)
/// assert!((energy - 4.935).abs() < 0.5);
/// ```
pub struct SchrodingerSolver<const N: usize> {
    /// Tableau du potentiel discret `V(xᵢ)` sur la grille spatiale.
    pub potential: [f32; N],
    /// Pas spatial `dx` (en mètres ou en unités réduites selon le système).
    pub dx: f32,
    /// Constante de Planck réduite `ħ` (en J·s ou unités réduites).
    pub hbar: f32,
    /// Masse de la particule `m` (en kg ou unités réduites).
    pub mass: f32,
}

impl<const N: usize> SchrodingerSolver<N> {
    /// Crée un nouveau solveur.
    ///
    /// # Arguments
    ///
    /// - `potential` — valeurs de `V(xᵢ)` aux `N` points de la grille.
    /// - `dx`        — pas spatial uniforme.
    /// - `hbar`      — constante de Planck réduite `ħ`.
    /// - `mass`      — masse de la particule.
    ///
    /// # Exemple
    ///
    /// ```rust
    /// use embedded_schrodinger::SchrodingerSolver;
    /// let solver = SchrodingerSolver::<64>::new([0.0f32; 64], 0.01, 1.0, 1.0);
    /// ```
    #[inline]
    pub fn new(potential: [f32; N], dx: f32, hbar: f32, mass: f32) -> Self {
        Self { potential, dx, hbar, mass }
    }

    /// Intègre la fonction d'onde de **gauche à droite** pour une énergie `energy` donnée.
    ///
    /// Applique la relation de récurrence issue de la discrétisation de l'équation
    /// de Schrödinger par **différences finies centrées** :
    ///
    /// ```text
    /// ψᵢ₊₁ = 2ψᵢ − ψᵢ₋₁ − (2m·dx²/ħ²)(E − V(xᵢ))ψᵢ
    /// ```
    ///
    /// Les conditions aux limites initiales sont celles d'un **puits infini** :
    /// `ψ(0) = 0`, `ψ(1) = 0.001` (impulsion de départ arbitraire).
    ///
    /// # Valeur de retour
    ///
    /// Valeur de `ψ` au dernier point de la grille (`ψ[N-1]`).
    /// Si la propagation diverge (NaN ou infini), la dernière valeur saine est renvoyée.
    pub fn integrate_wavefunction(&self, energy: f32, psi: &mut [f32; N]) -> f32 {
        psi[0] = 0.0;
        psi[1] = 0.001;

        // Facteur cinétique : (2·m·dx²) / ħ²
        let k_factor = (2.0 * self.mass * self.dx * self.dx) / (self.hbar * self.hbar);

        for i in 1..(N - 1) {
            psi[i + 1] = 2.0 * psi[i]
                - psi[i - 1]
                - k_factor * (energy - self.potential[i]) * psi[i];

            // Protection anti-divergence : évite de corrompre le calcul sur MCU
            if psi[i + 1].is_infinite() || psi[i + 1].is_nan() {
                return psi[i];
            }
        }

        psi[N - 1]
    }

    /// Normalise la fonction d'onde `ψ` selon la condition `∫|ψ|² dx = 1`.
    ///
    /// L'intégration numérique est effectuée par la **méthode des rectangles**.
    /// La racine carrée est calculée via
    /// [`embedded_f32_sqrt::sqrt`] (Newton-Raphson, aucun recours à `libm`).
    ///
    /// # Erreurs
    ///
    /// - `"Erreur mathématique lors du calcul de la norme (NaN/Infinity)"` si
    ///   la somme produit une valeur non-numérique.
    /// - `"Impossible de normaliser une fonction d'onde de norme nulle"` si la
    ///   norme calculée est nulle.
    pub fn normalize(&self, psi: &mut [f32; N]) -> Result<(), &'static str> {
        let mut sum = 0.0f32;

        for &val in psi.iter() {
            sum += val * val * self.dx;
        }

        let norm = sqrt(sum)
            .map_err(|_| "Erreur mathématique lors du calcul de la norme (NaN/Infinity)")?;

        if norm > 0.0 {
            for val in psi.iter_mut() {
                *val /= norm;
            }
            Ok(())
        } else {
            Err("Impossible de normaliser une fonction d'onde de norme nulle")
        }
    }

    /// Recherche un **état propre** (énergie + fonction d'onde) par la **méthode de tir**.
    ///
    /// L'algorithme effectue une dichotomie binaire sur l'intervalle `[e_min, e_max]`.
    /// À chaque itération, il intègre `ψ` et analyse la valeur frontière droite pour
    /// déterminer de quel côté se trouve la valeur propre recherchée.
    ///
    /// 30 itérations donnent une précision de l'ordre de `2⁻³⁰` de l'intervalle
    /// initial, soit la précision machine du type `f32`.
    ///
    /// # Arguments
    ///
    /// - `e_min` — borne basse de la plage de recherche (unité cohérente avec `hbar`/`mass`).
    /// - `e_max` — borne haute de la plage de recherche.
    /// - `psi`   — tampon mutable où sera écrite la fonction d'onde normalisée.
    ///
    /// # Valeur de retour
    ///
    /// - `Ok(energy)` — énergie propre convergée en unité choisie.
    /// - `Err(msg)`   — si la normalisation finale échoue.
    ///
    /// # Exemple
    ///
    /// ```rust
    /// use embedded_schrodinger::SchrodingerSolver;
    ///
    /// const N: usize = 200;
    /// let solver = SchrodingerSolver::new([0.0f32; N], 0.005, 1.0, 1.0);
    /// let mut psi = [0.0f32; N];
    ///
    /// let e1 = solver.find_eigenstate(1.0, 10.0, &mut psi).unwrap();
    /// // En unités atomiques, E₁ ≈ π²/2 ≈ 4.935
    /// assert!((e1 - 4.935).abs() < 0.3);
    /// ```
    pub fn find_eigenstate(
        &self,
        mut e_min: f32,
        mut e_max: f32,
        psi: &mut [f32; N],
    ) -> Result<f32, &'static str> {
        let mut e_mid = 0.0f32;

        for _ in 0..30 {
            e_mid = 0.5 * (e_min + e_max);
            let psi_boundary = self.integrate_wavefunction(e_mid, psi);

            // CORRECTION : logique de bissection corrigée.
            // ψ_boundary > 0 signifie qu'on est en-dessous de la valeur propre,
            // donc on relève e_min (et non e_max comme c'était le cas avant).
            if psi_boundary > 0.0 {
                e_min = e_mid;
            } else {
                e_max = e_mid;
            }
        }

        self.normalize(psi)?;

        Ok(e_mid)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests unitaires
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// En unités atomiques (ħ = 1, m = 1), dans un puits infini de largeur L,
    /// l'énergie du niveau fondamental est :
    ///   E₁ = π² / (2 L²)
    ///
    /// Pour L = 1.0 (N·dx = 1.0) : E₁ ≈ 4.9348
    #[test]
    fn test_energie_fondamentale_puits_infini() {
        const N: usize = 200;
        let dx = 1.0 / N as f32;
        let solver = SchrodingerSolver::new([0.0f32; N], dx, 1.0, 1.0);
        let mut psi = [0.0f32; N];

        let energy = solver
            .find_eigenstate(1.0, 15.0, &mut psi)
            .expect("find_eigenstate a échoué");

        // Référence analytique : π²/2 ≈ 4.9348
        let e_ref = core::f32::consts::PI * core::f32::consts::PI / 2.0;
        assert!(
            (energy - e_ref).abs() < 0.5,
            "Énergie fondamentale hors tolérance : {energy} vs {e_ref}"
        );
    }

    /// La norme de la fonction d'onde normalisée doit être égale à 1
    /// (à la précision numérique de l'intégration rectangulaire).
    #[test]
    fn test_normalisation_norme_unitaire() {
        const N: usize = 200;
        let dx = 1.0 / N as f32;
        let solver = SchrodingerSolver::new([0.0f32; N], dx, 1.0, 1.0);
        let mut psi = [0.0f32; N];

        solver
            .find_eigenstate(1.0, 15.0, &mut psi)
            .expect("find_eigenstate a échoué");

        let norm_sq: f32 = psi.iter().map(|&v| v * v * dx).sum();
        assert!(
            (norm_sq - 1.0).abs() < 1e-3,
            "La norme au carré devrait être ≈ 1, obtenu : {norm_sq}"
        );
    }

    /// Un potentiel constant non nul décale l'énergie propre.
    /// Pour V₀ = 1.0, E₁_décalée ≈ E₁_libre + V₀.
    #[test]
    fn test_potentiel_constant_decalage_energie() {
        const N: usize = 200;
        let dx = 1.0 / N as f32;

        let solver_libre = SchrodingerSolver::new([0.0f32; N], dx, 1.0, 1.0);
        let solver_offset = SchrodingerSolver::new([1.0f32; N], dx, 1.0, 1.0);

        let mut psi = [0.0f32; N];

        let e_libre = solver_libre
            .find_eigenstate(1.0, 15.0, &mut psi)
            .expect("Solver libre a échoué");

        let e_offset = solver_offset
            .find_eigenstate(2.0, 16.0, &mut psi)
            .expect("Solver décalé a échoué");

        assert!(
            (e_offset - e_libre - 1.0).abs() < 0.5,
            "Décalage en énergie incorrect : {e_offset} - {e_libre} ≠ 1.0"
        );
    }

    /// Vérifie que la normalisation renvoie une erreur sur un tableau nul.
    #[test]
    fn test_normalisation_vecteur_nul_retourne_erreur() {
        const N: usize = 50;
        let solver = SchrodingerSolver::new([0.0f32; N], 0.01, 1.0, 1.0);
        let mut psi = [0.0f32; N];

        let result = solver.normalize(&mut psi);
        assert!(result.is_err(), "Devrait échouer sur une fonction d'onde nulle");
    }

    /// Vérifie que `integrate_wavefunction` retourne une valeur finie
    /// pour une énergie physiquement raisonnable.
    #[test]
    fn test_integration_produit_valeur_finie() {
        const N: usize = 100;
        let dx = 0.01;
        let solver = SchrodingerSolver::new([0.0f32; N], dx, 1.0, 1.0);
        let mut psi = [0.0f32; N];

        let boundary = solver.integrate_wavefunction(5.0, &mut psi);
        assert!(
            boundary.is_finite(),
            "La valeur frontière doit être finie, obtenu : {boundary}"
        );
    }
}