# embedded-schrodinger
[![License: GPL-2.0](https://img.shields.io/badge/License-GPL--2.0-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![no_std](https://img.shields.io/badge/no__std-✓-green)](https://docs.rust-embedded.org/book/intro/no-std.html)
[![no libm](https://img.shields.io/badge/libm-✗-red)](https://crates.io/crates/libm)

Solveur de l'**équation de Schrödinger indépendante du temps en 1D** conçu
pour les systèmes **embarqués** (`no_std`, sans `libm`, sans liaison C,
sans `unsafe`).

---

## Méthode numérique

La résolution repose sur la **méthode de tir** (*shooting method*) :

1. Pour une énergie candidate `E`, la fonction d'onde `ψ` est propagée de
   gauche à droite via la relation de récurrence des **différences finies centrées** :
```text
   ψᵢ₊₁ = 2ψᵢ − ψᵢ₋₁ − (2m·dx²/ħ²)(E − V(xᵢ))ψᵢ
```
2. Une **dichotomie binaire** (30 itérations) ajuste `E` jusqu'à ce que
   `ψ(xₙ) ≈ 0` (condition aux limites du puits infini).
   La logique de signe est la suivante :
   - `ψ(xₙ) > 0` → `E` est trop faible → on relève `e_min`
   - `ψ(xₙ) < 0` → `E` est trop élevé → on abaisse `e_max`
3. La fonction d'onde convergée est **normalisée** : `∫|ψ|² dx = 1`.

---

## Caractéristiques

| Propriété         | Valeur                                         |
|-------------------|------------------------------------------------|
| `no_std`          | ✅ oui                                         |
| `libm`            | ❌ aucun recours                               |
| Liaison C         | ❌ aucune                                      |
| `unsafe`          | ❌ `#![forbid(unsafe_code)]`                   |
| Racine carrée     | Newton-Raphson via `embedded-f32-sqrt`         |
| Précision         | `f32` (IEEE 754 simple précision)              |
| Cibles validées   | Cortex-M33, Cortex-M4F et hôte                |

---

## Installation

Dans votre `Cargo.toml` :

```toml
[dependencies]
embedded-schrodinger = "0.1"
```

---

## Exemple rapide

```rust
use embedded_schrodinger::SchrodingerSolver;

const N: usize = 200;

// Puits de potentiel infini (V = 0 à l'intérieur)
let potential = [0.0f32; N];

// Unités atomiques : ħ = 1, m = 1, dx = L/N avec L = 1
let solver = SchrodingerSolver::new(potential, 1.0 / N as f32, 1.0, 1.0);

let mut psi = [0.0f32; N];
let energy = solver
    .find_eigenstate(1.0, 15.0, &mut psi)
    .expect("Aucune solution dans l'intervalle donné");

// E₁ = π²/2 ≈ 4.935 en unités atomiques pour L = 1
println!("Énergie fondamentale : {:.4}", energy);
```

---

## Dépendances

| Crate               | Rôle                                         |
|---------------------|----------------------------------------------|
| `embedded-f32-sqrt` | Racine carrée `f32` Newton-Raphson, `no_std` |

Aucune autre dépendance. Aucun accès au système d'exploitation.

---

## Licence

Ce logiciel est distribué sous les termes de la **GNU General Public License
version 2** (GPL-2.0-only).

Voir [https://www.gnu.org/licenses/old-licenses/gpl-2.0.html](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
pour le texte complet de la licence.

**Copyright (C) 2025  Jorge Andre Castro.**