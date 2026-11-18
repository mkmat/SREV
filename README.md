<pre>
    __   __   __     
   |__  |__| |__ \  /   (c) 2026 martin kroger
    __| |  \ |__  \/             marcelo carignano
</pre>

# SREV

self-returning excluded volume path

## start from a parameter file 

        perl SREV -dir=<directory-name> [-np=<cores>] [-sphere] [-silent] [-compile]

Options:

        -dir    : directory name containing SREV-parameters.txt (results are written into this directory) 
        -np     : number of cores used by mpirun (default: 1)
        -sphere : grows the SREV inside a sphere (default: inside a cubic box with periodic boundary conditions)
        -silent : do not generate log.lammps and further runtime information to stdout (default: verbose)
        -compile: recompile SREV. Only required if user-specific settings have been modified inside SREV

Examples: 

        perl SREV -dir=runs/pbc-smoke-001 
        perl SREV -dir=runs/pbc-smoke-001 -np=2 -sphere

## specify parameters on the command line 

        perl SREV -N=<beads> -alpha=<value> -bmax=<value> -phi=<value> [-seed=<integer>] [-np=<cores>] [-sphere] [-silent] [-compile]

Options: 

        -N      : number of beads (results are written into the calling directory)
        -alpha  : SRRW parameter alpha
        -bmax   : maximum SRRW bond length
        -phi    : target volume fraction
        -seed   : fixed seed value (default: random seed)
        -np     : number of cores used by mpirun (default: 1)
        -sphere : grows the SREV inside a sphere (default: inside a cubic box with periodic boundary conditions)
        -silent : do not generate log.lammps and further runtime information to stdout (default: verbose)
        -compile: recompile SREV. Only required if user-specific settings have been modified inside SREV

Examples: 

        perl SREV -N=5000 -alpha=1.1 -bmax=20.0 -phi=0.12 
        perl SREV -N=5000 -alpha=1.1 -bmax=20.0 -phi=0.12 -seed=12345 -np=10 -sphere
