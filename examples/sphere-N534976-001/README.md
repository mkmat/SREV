# Spherical N=534976 run 001

Full-length spherical-confinement test using the current source and the same
scientific parameters and fixed seed as the validated PBC calculation. 

Parameters:

- `N = 534976` bonds (`534977` beads)
- `alpha = 1.1`
- `bmax = 20.0`
- `phi = 0.12`
- `node_radius = 0.49`
- `fixed_seed = 12345`
- two MPI processes for LAMMPS

For this run only, the build copy selects
`L = L_for_sphere(); call SRRW_in_sphere` in place of the default PBC call.
The tracked source remains configured for PBC.

Results:

- spherical diameter `L = 161.291` (`L/2 = 80.6455`)
- wall time: `72.35 s`
- final minimum bead-bead distance: `0.991342`
- mean SRRW-to-SREV displacement: `1.030`
- measured SREV mean bond length: `2.375`

The maximum radius after final recentering is `103.8904` for the SRRW and
`103.9889` for the SREV. This is expected from the existing evolving-center-of-
mass definition of spherical confinement: it is the same caveat before and
after the reject-and-retry change.

OVITO files:

- `SREV-sphere-nonperiodic.dump` is the SREV dump with nonperiodic boundary
  flags for correct spherical coordination analysis.
- `SREV_sphere_CN_cutoff1.15.ovito` is the saved OVITO view, colored by
  coordination number with cutoff `1.15`, range `0` to `13`, periodic wrapping
  disabled, and the simulation-cell display hidden.
