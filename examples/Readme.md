Each directory in this folder originally contained only two files: 

1. SREV-parameters.txt (required by SREV)
2. README.md (for your information)

Upon calling 

    perl SREV -dir="examples/sphere-N534976-001" -sphere -np=8

the SREV code produced the following files: 

    backbone-unfolded.dump
    _bead_bead_contacts.txt
    _Cp_id.txt
    log.lammps
    _SREV_average_linker_bp.txt
    _SREV-bondlength.txt
    SREV.data
    SREV_sphere_CN_cutoff1.15.ovito
    SREV-unfolded.dump
    SREV-unfolded.xyz
    SRRW-unfolded.dump

Some of them are too large to be uploaded here. We include some of them. 

Upon calling 

    perl SREV -dir="examples/cube-N534976-001" -np=8

the SREV produced a similar set of files, also included here. Note that the two SREV-parameters.txt files in the two subdirectories are identical! The difference in the geometry is taken into account by the SREV argument -sphere. 

