# pp2dI8_Ber_Split

An extension of the FeatFlow application `pp2d`, a solver for the 2D incompressible
Navier–Stokes equations based on discrete projection and multigrid.

This version adds:

- higher-order **Bernstein mixed finite elements** (Q2/P1disc and Q3/P2disc),
- **operator splitting** (Glowinski type) within the fractional-step scheme,
- **stabilization** of the convection by jump terms or algebraic flux correction (AFC/FCT),
- **isoparametric elements** for curved boundaries.

All new code is in `src/`; the FEAT library is unchanged. The features are selected
through the usual parameter file `#data/pp2d.dat`, which keeps the original pp2d format.

## Build and run

The directory must sit inside a FeatFlow installation. Then:

```sh
make                      # build the executable
./pp2d-pc64-x86-linux     # run with #data/pp2d.dat and the mesh in #pre/
```
