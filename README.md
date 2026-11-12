# pp2dI8_Ber_Split

FeatFlow `pp2d` (2D incompressible Navier–Stokes, discrete projection, multigrid)
extended with **Bernstein Qk / P(k−1)disc mixed elements**, Glowinski operator
splitting, AFC/FCT stabilization of the convection and isoparametric boundary
elements. The FEAT library is untouched; all new code lives in `src/`.

## Features

| Feature | Selected by (`#data/pp2d.dat`) |
|---|---|
| Q1~/Q0 (original pp2d) | `IELT = 0..3` |
| Bernstein Q2 / P1disc | `IELT = 4` |
| Bernstein Q3 / P2disc | `IELT = 5` |
| One-step θ / fractional step (coupled) | `IFRSTP = 0 / 1` |
| Fractional step + Glowinski splitting (ISPLIT = 1, 2, 3) | `IFRSTP = 2, 3, 4` |
| Gradient-jump stabilization | `IJUMP = 2`, `DJUMP` (needs `IUPW = 0`) |
| AFC/FCT of the convection (Q1~, Q2, also with splitting) | `IAFC = 1` (Zalesak), `2` (test, α = 1), `CQAFC` |
| Isoparametric curved boundary elements (Q2, Q3) | `IISO = 1` |

`pp2d.dat` keeps the original pp2d structure. Additional options (`IAFC`, `CQAFC`,
`IISO`, …) are optional lines at the end of the file, read by name in any order;
unknown names produce a warning and are ignored. Old switches (`IBERN`, `IBNTST`,
`ISPLIT`) are still accepted.

Constraints checked in `GDAT` (init1.f):
- AFC needs `IJUMP = 0`, `IREFIN = 1`, `IPRECA = 1`, `IUPW = 0`, `UPSAM = 0`.
- `IISO = 1` needs `IJUMP = 0`.

## Element details

- Velocity: tensor-product Bernstein polynomials on the reference square, C⁰ by
  shared vertex/edge coefficients; global dofs `NVT + (k−1)·NMT + (k−1)²·NEL`.
- Pressure: mapped P(k−1) discontinuous, Bernstein basis on a reference triangle
  containing the square, `k(k+1)/2` dofs per element.
- Lumped Bernstein mass matrix has positive row sums for every k, so the
  lumped-mass projection `C = Bᵀ Ml⁻¹ B` of pp2d remains usable.
- Exact grid transfer for multigrid (`bnq3.f`, `bnpp.f`); maximum order `NNBK = 8`.

## Source overview (`src/`)

| File | Content |
|---|---|
| `bern.f` | Bernstein basis, dof numbering, cubature, own assembly |
| `bnpp.f` | pp2d integration: matrices, boundary conditions, MG transfer, forces |
| `bnq3.f` | Generic Qk element (EB20/EB30), Q3 dofs, generic MG transfer |
| `bniso.f` | Isoparametric boundary map (Gordon–Hall blending) |
| `bnafc.f`, `afcns.f` | AFC/FCT (scalar test and Navier–Stokes convection) |
| `bnstok.f` | Stand-alone Stokes solver / tests |
| `jump.f` | Jump stabilization |
| `init1.f` | Parameter input (`GDAT`, `GDKEY` for named optional lines) |
| others | original pp2d (`pp2d.f`, `nsdef.f`, `prostp.f`, `mgrout.f`, …) |

## Build

Must sit inside a FeatFlow tree (`FEATFLOW = ../../..` in the `Makefile`) with the
`feat2dI8` and `sysutils` libraries built:

```sh
make            # builds pp2d-$(ID), e.g. pp2d-pc64-x86-linux
```

## Run

From this directory, with the parameters in `#data/pp2d.dat` and the mesh in `#pre/`
(default: flow around a cylinder, `c2d`):

```sh
./pp2d-pc64-x86-linux
```

Output goes to `#data/` (protocol `PP2D.PRT`, solution `#DX*`), `#gmv/`, `#avs/`,
`#ns/` and `#points/`; these are not tracked in git.
