![GitHub Logo](logo.png)

# MagPot

**MagPot: A magnetoelastic potential generator for molecular dynamics simulations**

**Authors:** P.Nieves, I.Korniienko, S.Arapan, A.Fraile, R.Iglesias and D.Legut

MagPot is a Python program that generates Bethe–Slater parameters for magnetic interactions in spin–lattice molecular dynamics simulations using the **SPIN** package of **LAMMPS**.

**Current release:** [MagPot v1.0](https://github.com/pnieves2019/magpot/releases/tag/v1.0)  
**License:** [MIT](LICENSE)

## Overview

MagPot converts a compact set of structural, elastic, magnetic and magnetoelastic material properties into LAMMPS-ready pair coefficients for:

- the Heisenberg exchange interaction, described by the radial function \(J(r)\);
- the two-ion Néel interaction, described by the dipolar and quadrupolar radial functions \(l(r)\) and \(q(r)\).

A central feature of MagPot is its **modular architecture**. The program generates the magnetic contribution to the spin–lattice model; it does **not** generate the spin-independent reference potential describing the underlying atomic interactions. The reference potential must be supplied separately and combined with the MagPot-generated magnetic interactions in LAMMPS, for example using `pair_style hybrid/overlay`.

This approach allows the magnetic contribution to be used with different compatible reference-potential models, provided that the combined model is physically appropriate and does not double-count magnetic contributions.

## Supported crystal structures and magnetic orders

Version 1.0 supports the following cubic lattice and collinear magnetic-order combinations:

| Crystal structure | Magnetic order | Configuration |
|---|---|---|
| Simple cubic (SC) | FM | Ferromagnetic |
| Simple cubic (SC) | AFM | Antiferromagnetic |
| Body-centred cubic (BCC) | FM | Ferromagnetic |
| Body-centred cubic (BCC) | AFM | G-type antiferromagnetic |
| Face-centred cubic (FCC) | FM | Ferromagnetic |
| Face-centred cubic (FCC) | AFM | Type-II antiferromagnetic |

**Important note on FCC-AFM (type II).** For this configuration, MagPot retains the six antiparallel second-neighbour interactions, which form the nearest-neighbour shell of each of four interpenetrating simple-cubic networks. First-neighbour interactions and interactions between different networks are neglected. The corresponding LAMMPS structure must use four atom types to identify the four networks. These atom types identify crystallographic networks, not spin directions; spin vectors must be assigned independently. See the User Manual and the article for details.

Version 1.0 is limited to the structures and magnetic orders listed above. It does not generate potentials for non-cubic structures or general non-collinear magnetic configurations.

## Repository contents

| File or directory | Description |
|---|---|
| `magpot.py` | Main MagPot Python program. |
| `input.dat` | Example input file. Edit this file for the material and magnetic configuration of interest. |
| `requirements.txt` | Python package requirements. |
| `Examples/` | Example calculations and LAMMPS input files. |
| `User_Manual.pdf` | User documentation and additional instructions. |
| `LICENSE` | MIT license. |

## Requirements

- Python 3
- NumPy 1.18.4 or later
- LAMMPS with the SPIN package to run spin–lattice simulations using the generated coefficients

MagPot itself generates the parameters and LAMMPS commands; LAMMPS is required for the subsequent simulations.

## Installation

Clone the repository and install the Python dependency:

```bash
git clone --branch v1.0 --single-branch https://github.com/pnieves2019/magpot.git
cd magpot
python3 -m pip install -r requirements.txt
```

Alternatively, download the source archive from the [v1.0 release page](https://github.com/pnieves2019/magpot/releases/tag/v1.0).

## Quick start

1. Open `input.dat` and replace the example values with the material properties for your system. The input file must be in the directory from which you run `magpot.py`.
2. Run MagPot:

   ```bash
   python3 magpot.py > results.dat
   ```

3. Inspect `results.dat`. It contains the generated Bethe–Slater parameters and the corresponding LAMMPS `pair_coeff` commands. Use these magnetic coefficients in your LAMMPS input script together with a separately specified reference potential.

## Input parameters

After comment lines (beginning with `#`) are ignored, `input.dat` must contain **11 entries in the following order**:

| Entry | Meaning | Units / allowed values |
|---|---|---|
| `p` | Output format | `lmp` |
| `c` | Crystal structure | `sc`, `bcc` or `fcc` |
| `m` | Magnetic order | `fm` or `afm` |
| `rNN` | First-neighbour distance of the underlying magnetic cubic lattice at the equilibrium volume of the reference system | Å |
| `B0` | Bulk modulus of the spin-independent reference system | GPa |
| `Tord` | Curie temperature for FM order or Néel temperature for AFM order | K |
| `ws` | Spontaneous volume magnetostriction at zero temperature | Dimensionless |
| `K1` | First magnetocrystalline-anisotropy constant at zero temperature | J m⁻³ |
| `eta` | Relative pressure derivative of `K1`, \(\eta = K_1^{-1}\,\partial K_1/\partial P\) | GPa⁻¹ |
| `b1` | First anisotropic magnetoelastic constant | MPa |
| `b2` | Second anisotropic magnetoelastic constant | MPa |

Use the sign conventions and definitions described in the User Manual and the accompanying article. In particular, `rNN` and `B0` refer to the spin-independent reference system, while the remaining entries describe target magnetic and magnetoelastic properties. For FCC-AFM type II, MagPot internally converts the FCC first-neighbour distance to the second-neighbour interaction distance used by the model.

## Using the generated coefficients in LAMMPS

The generated exchange and Néel coefficients are intended to be combined additively with a separately supplied reference potential. A schematic pair-style declaration is:

```lammps
atom_style spin
pair_style hybrid/overlay <reference-pair-style> spin/exchange <Rc> spin/neel <Rc>
```

Replace `<reference-pair-style>` and `<Rc>` with the appropriate reference-potential style and the cutoff value printed by MagPot. The reference-potential `pair_coeff` commands must be provided independently, together with the magnetic `pair_coeff` commands generated by MagPot.

The SPIN pair styles require a spin-capable atom style (`atom_style spin`). For FCC-AFM type II, follow the four-atom-type convention described above and in the User Manual. Consult the LAMMPS documentation for the syntax and requirements of the reference potential and SPIN package.

## Examples and documentation

The [`Examples/`](Examples/) directory contains example calculations and LAMMPS input files. The [User Manual](User_Manual.pdf) describes the input format and program operation in more detail.

## Article and citation

If you use MagPot in published work, please cite the associated article/preprint and identify the software version used.

**Preprint submitted to arXiv**  
P. Nieves, I. Korniienko, S. Arapan, A. Fraile, R. Iglesias and D. Legut, *MagPot: A magnetoelastic potential generator for molecular dynamics simulations*.

**Software version**  
MagPot v1.0, [GitHub release](https://github.com/pnieves2019/magpot/releases/tag/v1.0).

## License

MagPot is distributed under the [MIT License](LICENSE). See the `LICENSE` file for the full license text.

