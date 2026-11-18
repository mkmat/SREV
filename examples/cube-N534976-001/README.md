# PBC N=534976 run 001

Full-length reproducible PBC test corresponding to the remembered configuration
name `SRRW-PBC-L=130-N=534976-phi=0.12-alpha=1.1-bmax=20-r=0.49-ru=10`.

Parameters read by the Fortran program:

- `N = 534976`
- `alpha = 1.1`
- `bmax = 20.0`
- `phi = 0.12`
- `fixed_seed = 12345`
- `bond_and_angle_forces = F`
- `first_config = T`
- `node_radius = 0.49` (defined in the source)
- MPI processes: 2

The source derives `L = 130` from `N`, `phi`, and `node_radius`.

The historical `ru=10` label is not represented by a parameter in this source
and is therefore not applied in this run.

## Result

Completed successfully on 2026-11-17 using two MPI processes.

- Derived box length: `L = 130.000`
- SREV coordinates: 534,977
- Final reported minimum separation: 0.996387
- Mean SRRW-to-SREV displacement: 0.876
- Measured SREV mean bond length: 2.334
- Contact pairs: 19,169,190
- Wall time: 38.62 seconds

Periodic minimum-image validation found four bead pairs with center-to-center
distance strictly below 1.0. They are listed in `pairs-distance-lt-1.txt`.

The OVITO session `SREV_CN_cutoff1.15.ovito` displays the configuration wrapped
into the periodic box and colored by coordination number using cutoff 1.15. The
blue-white-red color scale spans coordination numbers 0 through 13.

## Validation after disabling TAD code

The updated single-type program was run again on 2026-11-17 with the same
parameters, fixed seed, and two MPI processes. It completed successfully in
34.39 seconds. `SREV-unfolded.xyz` and `_bead_bead_contacts.txt` have exactly
the same SHA-256 checksums as the original validated run, confirming that the
coordinates and contact results are unchanged. `SREV.data` intentionally now
contains one generic atom type.

## Validation after removing beta

The independent `beta` parameter was removed because the validated model uses
the same exponent for both distributions. Repeating this full run with
`alpha=1.1` completed successfully in 33.83 seconds. The coordinate and contact
list SHA-256 checksums are unchanged.

