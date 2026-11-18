! (c) 2026 SREV is copyrighted by the copyright holders Martin Kröger and Marcelo Carignano. 
! It may be freely used for educational and research purposes by non-profit institutions and US government agencies only. 
! All other organizations may use SREV for evaluation purposes only, and any further uses will require prior written approval. 
! This software may not be sold or redistributed without prior written approval. 
! Copies of the software may be made by a user provided that copies are not sold or distributed, 
! and provided that copies are used under the same terms and conditions as agreed to in this paragraph.
! As research software, this code is provided on an "as is'' basis without warranty of any kind, either expressed or implied. 
! The downloading, or executing any part of this software constitutes an implicit agreement to these terms. 
! These terms and conditions are subject to change at any time without prior notice.

! use the SREV script to compile and run this code.

! compile: gfortran -Ofast -ffree-line-length-none SREV-19nov26.f90 -o SREV.ex
! or     : gfortran -O -ffree-line-length-none code.f90 -o code.ex
!
! LAMMPS installation note:
! The MPI launcher options and LAMMPS executable path used are specific to the authors' workstation.
! Users must adapt them inside SREV to their MPI and LAMMPS installations.

! 2025 ETH Zurich. martin kroger, mk@mat.ethz.ch (do not re-distribute)

module shared
  implicit none

  ! SRRW parameters
  integer              :: N                                 ! number of beads
  real*8               :: alpha                             ! SRRW parameter
  real*8               :: bmax                              ! SRRW parameter
  real*8               :: phi                               ! volume fraction (L is calculated from phi below)  
  integer              :: sphere                            ! 1 for sphere, 0 for rectangular box with periodic boundary conditions
  real*8               :: node_radius = 0.49                ! bead radius used for volume fraction and excluded volume
  integer              :: fixed_seed = 0                    ! if fixed_seed > 0, the value is used, otherwise it is random  16 jan 2025
  logical              :: save_contact_map = .true.         ! save contact map 
  logical              :: first_config = .true.             ! added 30 jan 2025 (if .true. id-bp conversion factor is saved, if .false. it is used)
  logical              :: silent

  real*8, dimension(:,:), allocatable :: xu                 ! unfolded coordinates of SRRW or SREV (if generated)
  integer, dimension(:), allocatable  :: type               ! bead type used in dump files; zero by default
  real*8               :: L,chi,mean_bondlen_for_this_SRRW_model,conversion_bondlen
  integer              :: seed_size
  integer, allocatable :: seed(:)
  character(len=255)   :: launcher
    
  ! postprocessing parameters
  real*8               :: max_distance  = 5.0
  integer              :: distance_bins = 300

  ! N: number of bonds
  ! L: box size
  ! xu: unfolded coordinates xu(0:N,:)
  ! bp_for_id(id) creates bp number for bead id in {0,1,..,N}
  ! id_for_bp(bp) creates bead id for bp number

  ! available locally, inside some subroutines:
  ! bb: backbone vectors bb(1:bbid,:)
  ! b:  bond vectors b(1:N,:)
  ! xubb: unfolded backbone coordinates xubb(0:bbid,:)

contains

  ! ================================================================
  subroutine read_parameters
    implicit none
    logical  :: EX
    inquire(file='SREV-parameters.txt',exist=EX) 
    if (.not.EX) stop 'missing SREV-parameters.txt'
    open(2,file='SREV-parameters.txt'); rewind(2)
    read(2,*) N,alpha,bmax,phi,fixed_seed,silent,first_config
    inquire(file='.silent',exist=EX); if (EX) silent=.true.
    inquire(file='.sphere',exist=EX); if (EX) sphere=1
    close(2)
    return
  end subroutine read_parameters

  ! ================================================================
  subroutine read_launcher
    open(2,file='.launcher'); read(2,'(A)') launcher; close(2); launcher=trim(launcher)
  return
  end

  ! ================================================================
  subroutine setup_seed_and_constants
    implicit none
    ! allocate xu and type
    call read_parameters
    call read_launcher
    allocate(xu(0:N,3))
    allocate(type(0:N)); type = 0
    ! create random seed array, use fixed_seed if larger zero, and keep seed value for later use 
    call random_seed(size=seed_size); allocate(seed(seed_size));
    if (fixed_seed.gt.0) then
       seed = fixed_seed
       call random_seed(put=seed)
    endif
    call random_seed(get=seed)
    chi = 1.D0-bmax**(-1.D0-alpha)  ! just an abbreviation
    return
  end subroutine setup_seed_and_constants

  ! ===================================================================
  function vectorlength(vector)
    implicit none
    real*8   :: vector(3),vectorlength
    vectorlength = dsqrt(dot_product(vector,vector))
    return
  end function vectorlength

  ! ===================================================================
  function random_unit_vector()
    implicit none
    real*8   :: random_unit_vector(3),z,phi,sz
    real     :: myrand(2)
    real*8   :: pi=dacos(-1.D0)
    call random_number(myrand)
    z     = -1.D0+2.D0*myrand(1)
    phi   = 2.D0*pi*myrand(2)
    sz    = dsqrt(1.D0-z**2)
    random_unit_vector = [cos(phi)*sz,sin(phi)*sz,z]
    return
  end function random_unit_vector

  ! ===================================================================
  function new_bondlen()
    implicit none
    real*8   :: new_bondlen
    real     :: myrand
    call random_number(myrand)
    new_bondlen = dexp(-dlog(1.D0-chi*myrand)/(1.D0+alpha))
    return
  end function new_bondlen

  ! ===================================================================
  function new_bond_vector()
    implicit none
    real*8   :: new_bond_vector(3)
    new_bond_vector = new_bondlen()*random_unit_vector()
    return
  end function new_bond_vector

  ! ===================================================================
  function bp_for_id(id)
    integer           :: id
    integer(kind=4)   :: bp_for_id
    bp_for_id = int(id*(147 + (conversion_bondlen-0.98)/0.034))    
    return
  end function bp_for_id

  ! ===================================================================
  function id_for_bp(bp)
    integer           :: id_for_bp
    integer(kind=4)   :: bp
    id_for_bp = int(bp/dble(147 + (conversion_bondlen-0.98)/0.034))
    return
  end function id_for_bp

  ! ===================================================================
  subroutine SRRW_in_periodic_box
    implicit none
    real*8,  dimension(:,:), allocatable   :: b,bb
    real*8,  dimension(:,:), allocatable   :: xubb
    integer, dimension(:),   allocatable   :: typebb
    real*8   :: return_prob
    integer  :: id,bbid
    real     :: myrand
    ! if (.not.silent) print '(A)','subroutine SRRW_in_periodic_box ..'
    allocate(b(N,3),bb(N,3))
    xu    = 0.D0
    bbid  = 0
    id    = 0
    do while (id.lt.N)
       id = id + 1
       if (bbid.le.1) then
          return_prob = 0.D0
       else
          return_prob = vectorlength(bb(bbid,:))**(-alpha)/alpha
       endif
       call random_number(myrand)
       if (myrand.lt.return_prob) then
          b(id,:)     = -bb(bbid,:)
          bbid        = bbid - 1
       else
          bbid        = bbid + 1
          bb(bbid,:) = new_bond_vector()
          b(id,:)     = bb(bbid,:)
       endif
       xu(id,:) = xu(id-1,:) + b(id,:)
    enddo
    ! create and save backbone coordinates
    allocate(xubb(0:bbid,3))
    allocate(typebb(0:bbid))
    typebb = 1
    xubb(0,:) = xu(0,:)
    do id=1,bbid
       xubb(id,:) = xubb(id-1,:) + bb(id,:)
    enddo
    call dump_unfolded_walk(bbid,xubb,typebb,'backbone-unfolded.dump')
    return
  end subroutine SRRW_in_periodic_box

  subroutine SRRW_in_sphere
    implicit none
    real*8,  dimension(:,:), allocatable   :: b,bb
    real*8   :: com(3),com_try(3),b_try(3),xu_try(3)
    real*8   :: return_prob
    real*8,  dimension(:,:), allocatable   :: xubb
    integer, dimension(:),   allocatable   :: typebb
    integer  :: id,bbid,bbid_before
    real     :: myrand
    ! if (.not.silent) print '(A)','subroutine SRRW_in_sphere ..'
    allocate(b(N,3),bb(N,3))
    xu    = 0.D0
    bbid = 0
    com  = 0.D0
    id   = 0
    do while (id.lt.N)
       id = id + 1
       ! Spherical confinement by rejection sampling (not reflection):
       ! propose a complete SRRW move, including an ordinary return move when
       ! selected, and test its endpoint against the prospective center of
       ! mass. If the endpoint lies outside L/2, restore the backbone cursor
       ! and retry this same bead with a new random move. Thus an outside bead
       ! is never committed to xu, unlike the previous next-step backtracking
       ! treatment. This is the structured-loop equivalent of ERIC's retry.
       ! Note that the test uses the center of mass at the current growth step;
       ! later center-of-mass drift can move earlier beads beyond the final L/2.
       trial_step: do
          bbid_before = bbid

          if (bbid.le.1) then
             return_prob = 0.D0
          else
             return_prob = vectorlength(bb(bbid,:))**(-alpha)/alpha
          endif

          call random_number(myrand)
          if (myrand.lt.return_prob) then
             b_try = -bb(bbid,:)
             bbid  = bbid - 1
          else
             bbid       = bbid + 1
             bb(bbid,:) = new_bond_vector()
             b_try      = bb(bbid,:)
          endif

          xu_try  = xu(id-1,:) + b_try
          com_try = com + xu_try
          if (vectorlength(xu_try-com_try/dble(id+1)).gt.L/2.D0) then
             bbid = bbid_before
             cycle trial_step
          endif

          b(id,:)  = b_try
          xu(id,:) = xu_try
          com      = com_try
          exit trial_step
       enddo trial_step
    enddo
    ! enforce zero center of mass
    com = com/dble(N+1)
    do id=0,N
       xu(id,:)=xu(id,:)-com
    enddo
    ! create and save backbone coordinates
    allocate(xubb(0:bbid,3))
    allocate(typebb(0:bbid))
    typebb = 1
    xubb(0,:) = xu(0,:)
    do id=1,bbid
       xubb(id,:) = xubb(id-1,:) + bb(id,:)
    enddo
    call dump_unfolded_walk(bbid,xubb,typebb,'backbone-unfolded.dump')
    return
  end subroutine SRRW_in_sphere

  ! ===================================================================
  subroutine dump_unfolded_walk(local_N,local_xu,local_type,filename)
    implicit none
    character(len=*) :: filename
    integer  :: local_N,local_type(0:local_N),id
    real*8   :: local_xu(0:local_N,3)
    ! if (.not.silent) print '(A)','subroutine dump_unfolded_walk '//trim(filename)//' ..'
    open(2,file=trim(filename)); rewind(2)
    write(2,'(A,/,I1)') 'ITEM: TIMESTEP',0
    write(2,'(A,/,I8)') 'ITEM: NUMBER OF ATOMS',local_N+1
    write(2,'(A)') 'ITEM: BOX BOUNDS pp pp pp'
    write(2,*) -L/2,L/2
    write(2,*) -L/2,L/2
    write(2,*) -L/2,L/2
    write(2,'(A)') 'ITEM: ATOMS id type xu yu zu'
    do id=0,local_N
       write(2,'(I8,1x,I1,1x,3(F12.6,1x))') id+1,local_type(id),local_xu(id,:)    ! id+1 as lammps does not accept id 0
    enddo
    close(2)
    return
  end subroutine dump_unfolded_walk

  ! ===================================================================
  subroutine dump_folded_walk(local_N,local_xu,local_type,filename)
    implicit none
    character(len=*) :: filename
    integer  :: local_N,local_type(0:local_N),id
    real*8   :: local_xu(0:local_N,3)
    ! if (.not.silent) print '(A)','subroutine dump_folded_walk '//trim(filename)//' ..'
    open(2,file=trim(filename)); rewind(2)
    write(2,'(A,/,I1)') 'ITEM: TIMESTEP',0
    write(2,'(A,/,I8)') 'ITEM: NUMBER OF ATOMS',local_N+1
    write(2,'(A)') 'ITEM: BOX BOUNDS pp pp pp'
    write(2,*) -L/2,L/2
    write(2,*) -L/2,L/2
    write(2,*) -L/2,L/2
    write(2,'(A)') 'ITEM: ATOMS id type xu yu zu'
    do id=0,local_N
       write(2,'(I8,1x,I1,1x,3(F12.6,1x))') id+1,local_type(id),local_xu(id,:)-L*anint(local_xu(id,:)/L)
    enddo
    close(2)
    return
  end subroutine dump_folded_walk

  ! =================================================================== 
  subroutine save_mean_bond_length
    implicit none
    integer           :: bps_for_node = 147
    integer           :: id
    real*8            :: bondlen,mean_bondlen
    integer(kind=4)   :: bpsum
    mean_bondlen = 0.D0
    bpsum        = 0
    do id=1,N
       bondlen      = vectorlength(xu(id,:)-xu(id-1,:))
       mean_bondlen = mean_bondlen + bondlen
       bpsum        = bpsum + int((bondlen-0.98)/0.034)
    enddo
    mean_bondlen = mean_bondlen/dble(N)
    if (.not.silent) print '(A40,F10.3)','SRRW mean bond length (ensemble)',mean_bondlen_for_this_SRRW_model
    if (.not.silent) print '(A40,I10)',  'SRRW bps per bond',bp_for_id(1)-147
    if (.not.silent) print '(A40,I10)',  'SRRW bps per node+bond',bp_for_id(1)
    if (.not.silent) print '(A50,F10.3)','id-bp conversion bond length used',conversion_bondlen  
    if (.not.silent) print '(A40,F10.3)','SREV mean bond length (measured)',mean_bondlen
    if (.not.silent) print '(A40,I10)',  'POST-SREV bps per bond',bpsum/N
    if (.not.silent) print '(A40,I10)',  'POST-SREV bps per node+bond',bpsum/N + bps_for_node
    open(13,file='_SREV_average_linker_bp.txt')
    write(13,*) bpsum/dble(N)    ! N bonds, N+1 nodes
    close(13)

    if (first_config) then                         
       conversion_bondlen = mean_bondlen
       open(22,file='_SREV-bondlength.txt'); rewind(22)
       write(22,*) conversion_bondlen
       write(22,*) N,alpha,bmax,phi
       close(22)
       ! if (.not.silent) print '(A)','## you may next choose first_config = .false.'
    endif

    return
  end subroutine save_mean_bond_length

  ! =================================================================== 
  subroutine create_contacts_brute_force
    implicit none
    real*8   :: contact_distance = 3.5
    integer  :: bead1,bead2
    real*8   :: distance
    if (.not.save_contact_map) return
    open(88,file='_bead_bead_contacts_brute_force.txt'); rewind(88)
    do bead1=0,N
       do bead2=bead1+1,N
          distance = vectorlength(xu(bead1,:)-xu(bead2,:))
          if (distance.le.contact_distance) write(88,*) min(bead1,bead2),max(bead1,bead2)
       enddo
    enddo
    close(88)
    return
  end subroutine create_contacts_brute_force

  ! =================================================================== 
  subroutine save_Cp_id_fast
    implicit none
    real*8            :: contact_distance = 3.5
    integer           :: max_bin
    integer           :: bin
    real*8            :: distance,Cp_value,lo(3),hi(3)
    integer(kind=4), dimension(:), allocatable   :: Cp,count_Cp
    ! neighborlist
    integer  :: cells(3),cellno(3),id,bead1,bead2
    real*8   :: cellsize(3)
    integer, dimension(:), allocatable        :: nextbead
    integer, dimension(:,:,:), allocatable    :: firstbead
    integer  :: ix,iy,iz
    integer  :: jx,jy,jz
    ! if (.not.silent) print '(A)','subroutine save_Cp_id_fast ..'
    max_bin = N/2
    if (.not.silent) print '(A30,I8)','Cp bins ',max_bin
    allocate(Cp(0:max_bin)); Cp = 0

    ! I first create neighbor lists using cells of size > contact_distance
    lo = +1.D30
    hi = -1.D30
    do id=0,N
       lo = min(lo,xu(id,:))
       hi = max(hi,xu(id,:))
    enddo
    if (.not.silent) print '(A30,3(F10.3,1x))','lo ',lo
    if (.not.silent) print '(A30,3(F10.3,1x))','hi ',hi
    ! add artificial empty shell to simplify algorithm below
    lo = lo - 2*contact_distance
    hi = hi + 2*contact_distance
    cells = int((hi-lo)/contact_distance)
    cells = min(200,cells)
    cellsize = (hi-lo)/cells
    if (.not.silent) print '(A30,3(I10,1x))','cells ',cells
    if (.not.silent) print '(A30,3(F10.3,1x))','cellsize ',cellsize
    allocate(nextbead(0:N))
    allocate(firstbead(cells(1),cells(2),cells(3)))
    nextbead = -1
    firstbead = -1
    do id=0,N
       cellno = 1+int((xu(id,:)-lo)/cellsize)
       nextbead(id) = firstbead(cellno(1),cellno(2),cellno(3))
       firstbead(cellno(1),cellno(2),cellno(3)) = id
    enddo

    if (save_contact_map) open(88,file='_bead_bead_contacts.txt'); rewind(88)

    ! now scan through neighbor list and fill the Cp field
    do ix=2,cells(1)-1
       do iy=2,cells(2)-1   ! ranges due to artificial shell
          do iz=2,cells(3)-1
             bead1 = firstbead(ix,iy,iz)
             do while (bead1.ne.-1)
                ! intra-cell
                bead2 = nextbead(bead1)
                do while (bead2.ne.-1)
                   ! (bead1,bead2) pair available
                   distance = vectorlength(xu(bead1,:)-xu(bead2,:)) ! bead-bead distance
                   bin = abs(bead2-bead1)
                   if (distance.le.contact_distance) then
                      if (bin.le.max_bin) Cp(bin) = Cp(bin) + 1
                      if (save_contact_map) write(88,*) min(bead1,bead2),max(bead1,bead2)
                   endif
                   bead2 = nextbead(bead2)
                enddo
                ! inter-cell
                do jx= 0,1
                   do jy=-1,1
                      do jz=-1,1
                         if (jx.eq.1.or.jx.eq.0.and.jy.eq.1.or.jx.eq.0.and.jy.eq.0.and.jz.eq.1) then   ! picks 13 relevant cases
                            bead2 = firstbead(ix+jx,iy+jy,iz+jz)
                            do while (bead2.ne.-1)
                               ! (bead1,bead2) pair available
                               distance = vectorlength(xu(bead1,:)-xu(bead2,:)) ! bead-bead distance
                               bin = abs(bead2-bead1)
                               if (distance.le.contact_distance) then
                                  if (bin.le.max_bin) Cp(bin) = Cp(bin) + 1
                                  if (save_contact_map) write(88,*) min(bead1,bead2),max(bead1,bead2)
                               endif
                               bead2 = nextbead(bead2)
                            enddo
                         endif
                      enddo
                   enddo
                enddo
                bead1 = nextbead(bead1)
             enddo
          enddo
       enddo
    enddo

    if (save_contact_map) close(88)
    if (save_contact_map) then
        if (.not.silent) print '(A)','created _bead_bead_contacts.txt'
    endif

    ! obtain the Cp normalization factors: count_Cp
    if (.not.silent) print '(A)','obtain normalization factors'
    allocate(count_Cp(0:max_bin)); count_Cp = 0
    do bin = 1,max_bin
       count_Cp(bin) = count_Cp(bin) + N+1 - bin
    enddo

    open(2,file='_Cp_id.txt'); rewind(2)
    do bin=0,max_bin
       if (count_Cp(bin).gt.0) then
          Cp_value = dble(Cp(bin))/dble(count_Cp(bin))
          write(2,*) bin,Cp_value
          ! write(42,*) bin,Cp_value  ! marcelo
       endif
    enddo
    return
  end subroutine save_Cp_id_fast

  ! =================================================================== 
  subroutine save_Cp_id_logarithmic
    implicit none
    real*8            :: contact_distance = 3.5
    real*8            :: base = 1.05            ! the more close the base is to unity, the more bins are generated
    integer           :: max_bin
    integer           :: bin
    real*8            :: distance,Cp_value,lo(3),hi(3)
    integer(kind=4), dimension(:), allocatable   :: Cp,count_Cp
    integer(kind=4)   :: id_distance
    ! neighborlist
    integer  :: cells(3),cellno(3),id,bead1,bead2
    real*8   :: cellsize(3)
    integer, dimension(:), allocatable        :: nextbead
    integer, dimension(:,:,:), allocatable    :: firstbead
    integer  :: ix,iy,iz
    integer  :: jx,jy,jz
    ! if (.not.silent) print '(A)','subroutine save_Cp_id_logarithmic ..'
    max_bin = nint(dlog(dble(N))/dlog(base))                ! bin = log(id_distance)/log(base)
    if (.not.silent) print '(A30,I8)','Cp bins ',max_bin
    allocate(Cp(0:max_bin)); Cp = 0

    ! I first create neighbor lists using cells of size > contact_distance
    lo = +1.D30
    hi = -1.D30
    do id=0,N
       lo = min(lo,xu(id,:))
       hi = max(hi,xu(id,:))
    enddo
    if (.not.silent) print '(A30,3(F10.3,1x))','lo ',lo
    if (.not.silent) print '(A30,3(F10.3,1x))','hi ',hi
    ! add artificial empty shell to simplify algorithm below
    lo = lo - 2*contact_distance
    hi = hi + 2*contact_distance
    cells = int((hi-lo)/contact_distance)
    cells = min(200,cells)
    cellsize = (hi-lo)/cells
    if (.not.silent) print '(A30,3(I10,1x))','cells ',cells
    if (.not.silent) print '(A30,3(F10.3,1x))','cellsize ',cellsize
    allocate(nextbead(0:N))
    allocate(firstbead(cells(1),cells(2),cells(3)))
    nextbead = -1
    firstbead = -1
    do id=0,N
       cellno = 1+int((xu(id,:)-lo)/cellsize)
       nextbead(id) = firstbead(cellno(1),cellno(2),cellno(3))
       firstbead(cellno(1),cellno(2),cellno(3)) = id
    enddo

    ! now scan through neighbor list and fill the Cp field
    do ix=2,cells(1)-1
       do iy=2,cells(2)-1   ! ranges due to artificial shell
          do iz=2,cells(3)-1
             bead1 = firstbead(ix,iy,iz)
             do while (bead1.ne.-1)
                ! intra-cell
                bead2 = nextbead(bead1)
                do while (bead2.ne.-1)
                   ! (bead1,bead2) pair available
                   distance = vectorlength(xu(bead1,:)-xu(bead2,:)) ! bead-bead distance
                   id_distance = abs(bead2-bead1)
                   bin = nint(dlog(dble(id_distance))/dlog(base))
                   if (distance.le.contact_distance) Cp(bin) = Cp(bin) + 1
                   bead2 = nextbead(bead2)
                enddo
                ! inter-cell
                do jx= 0,1
                   do jy=-1,1
                      do jz=-1,1
                         if (jx.eq.1.or.jx.eq.0.and.jy.eq.1.or.jx.eq.0.and.jy.eq.0.and.jz.eq.1) then   ! picks 13 relevant cases
                            bead2 = firstbead(ix+jx,iy+jy,iz+jz)
                            do while (bead2.ne.-1)
                               ! (bead1,bead2) pair available
                               distance = vectorlength(xu(bead1,:)-xu(bead2,:)) ! bead-bead distance
                               id_distance = abs(bead2-bead1)
                               bin = nint(dlog(dble(id_distance))/dlog(base))
                               if (distance.le.contact_distance) Cp(bin) = Cp(bin) + 1
                               bead2 = nextbead(bead2)
                            enddo
                         endif
                      enddo
                   enddo
                enddo
                bead1 = nextbead(bead1)
             enddo
          enddo
       enddo
    enddo

    ! obtain the Cp normalization factors: count_Cp
    if (.not.silent) print '(A30)','obtain normalization factors ..'
    allocate(count_Cp(0:max_bin)); count_Cp = 0
    do id_distance = 1,N
       bin = nint(dlog(dble(id_distance))/dlog(base))
       count_Cp(bin) = count_Cp(bin) + N+1 - id_distance   ! N contributions for id_distance = 1, N-1 for id_distance = 2 and so on
    enddo

    open(2,file='_Cp_id_logarithmic.txt'); rewind(2)
    do bin=0,max_bin
       if (count_Cp(bin).gt.0) then
          id_distance = nint(base**bin)
          Cp_value = dble(Cp(bin))/dble(count_Cp(bin))
          write(2,*) id_distance,Cp_value
          write(43,*) id_distance,Cp_value  ! marcelo
       endif
    enddo
    return
  end subroutine save_Cp_id_logarithmic

  subroutine create_SREV_via_lammps
    implicit none
    real*8   :: small = 0.05         ! keep
    real*8   :: soft = 15.0          ! keep
    integer  :: id,atomtypes
    logical  :: EX
    real*8   :: mindistance
    real*8   :: mean_SRRW_SREV_displacement,dx(3)
    logical  :: lammps_not_converged
    real*8, dimension(:,:), allocatable    :: xu_init
    character(len=255) :: silent_string
    allocate(xu_init(0:N,3))
    xu_init   = xu
    atomtypes = 1
    ! initial small displacement
    do id=1,N; xu(id,:)=xu(id,:)+small*random_unit_vector(); enddo

    ! create data file
    open(2,file='config.data'); rewind(2)
    write(2,'(A)') 'SREV 2026'
    write(2,'(A)') ' '
    write(2,'(I10,A)') N+1,' atoms'
    write(2,'(I10,A)') atomtypes,' atom types'
    write(2,'(A)') ' '
    write(2,'(2(F20.8,1x),A)') -L/2.D0,L/2.D0,' xlo xhi'
    write(2,'(2(F20.8,1x),A)') -L/2.D0,L/2.D0,' ylo yhi'
    write(2,'(2(F20.8,1x),A)') -L/2.D0,L/2.D0,' zlo zhi'
    write(2,'(A)') ' '
    write(2,'(A)') 'Atoms # atomic'
    write(2,'(A)') ' '
    do id=0,N
       write(2,'(I10,1x,I6,1x,3(F15.8,1x),A)') 1+id,1,xu(id,:),' 0 0 0'
    enddo
    close(2)

    if (silent) then
        silent_string = '';
    else
        silent_string = 'print "# minimum distance $r (amount $g at first bin, $counts counts, $problems problems)\n";'
    endif

    ! create temporary perl script
    open(2,file='_analyze_gr.pl'); rewind(2)
    write(2,'(10(A,/))') '#! /usr/bin/perl','open(G,"<gr.txt"); $line=<G>; $line=<G>; $line=<G>; $line=<G>;', &
         '$r=0; $problems = 0; $c=0; $firstbin=0; $secondbin=0;', &
         'while (!eof(G)) { $c+=1; $line=<G>; chomp $line; ($i,$r,$g,$counts)=split(/ /,$line);', &
         'if (($c eq 1)&($g>10000)) { $firstbin=1; }; if (($c>1)&($c<10)&($g>0)) { $secondbin=1; }; $problems += $g; }', &
         '`rm -f _stop`; if (($firstbin eq 1)&($secondbin eq 0)) { open(M,">_stop"); print M "1\n"; close(M); exit; };', &
         'open(G,"<gr.txt"); $line=<G>; $line=<G>; $line=<G>; $line=<G>;', &
         'while (!eof(G)) { $line=<G>; chomp $line; ($i,$r,$g,$counts)=split(/ /,$line); if ($g>0) {', &
         silent_string, & 
         'open(M,">_mindistance"); print M "$r\n"; close(M); exit; }; };'
    close(2)

    ! create lammps in-file
    open(2,file='_lammps.in'); rewind(2)
    write(2,'(A,F10.4)') 'variable PotEnergy equal pe'
    write(2,'(A)') 'units lj'
    write(2,'(A)') 'atom_style atomic'
    write(2,'(A)') 'read_data config.data'
    write(2,'(A)') 'neigh_modify delay 0 check yes one 100000 page 1000000'
    write(2,'(A)') 'mass * 1'
    write(2,'(A)') 'timestep 0.005'
    write(2,'(A,F10.3)') 'pair_style soft ',2.0*node_radius*1.03
    write(2,'(A,F10.3)') 'pair_coeff * * ',soft
    write(2,'(A)') 'fix 2 all nvt temp 0.001 0.001 1.0'
    write(2,'(A)') 'thermo 500'
    write(2,'(A)') 'thermo_style custom time pe v_PotEnergy step temp'
    write(2,'(A)') 'compute gr all rdf 200 1 1'
    write(2,'(A)') 'fix 1 all ave/time 1 1 1 c_gr[*] file gr-init.txt mode vector'
    write(2,'(A)') 'run 0'
    write(2,'(A)') '# variable mystep equal step'        
    write(2,'(A)') '# variable mytemp equal step'        
    write(2,'(A)') '# fix myprint all print 1 "${mystep} ${mytemp}" file _step_temperature.txt screen no'   
    write(2,'(A)') 'unfix 1'
    write(2,'(A)') 'uncompute gr'
    write(2,'(A)') 'shell rm -f SREV-unfolded.dump _mindistance gr-init.txt'
    write(2,'(A)') 'compute pairdist all pair/local dist'
    write(2,'(A)') 'compute mindist all reduce min c_pairdist inputs local'
    write(2,'(A)') 'variable mymindist equal c_mindist'
    write(2,'(A)') 'thermo_style custom time pe v_PotEnergy v_mymindist step temp'
    write(2,'(A)') 'run 0'
    write(2,'(A)') ' '
    write(2,'(A)') 'variable stuck equal 200'
    write(2,'(A)') 'label loop'
    write(2,'(A)') 'variable a loop ${stuck}'
    write(2,'(A)') 'variable LastPotEnergy equal ${PotEnergy}'
    write(2,'(A)') 'run 200'
    write(2,'(A)') 'compute gr all rdf 200 1 1'
    write(2,'(A)') 'fix 1 all ave/time 1 1 1 c_gr[*] file gr.txt mode vector'
    write(2,'(A)') 'run 0'
    write(2,'(A)') 'shell perl ./_analyze_gr.pl'
    write(2,'(A)') 'shell rm -f _stop'
    write(2,'(A)') 'unfix 1'
    write(2,'(A)') 'uncompute gr'
    write(2,'(A)') 'variable DiffEnergy equal abs(${PotEnergy}-${LastPotEnergy})'
    write(2,'(A)') 'if "${a}>=${stuck}" then "jump _lammps.in finish"'
    write(2,'(A,F7.4,A)') 'if "${mymindist}>',2*node_radius,'" then "jump _lammps.in break"'
    write(2,'(A)') 'next a'
    write(2,'(A)') 'jump _lammps.in loop'
    write(2,'(A)') 'label break'
    write(2,'(A)') 'variable a delete'
    write(2,'(A)') ' '
    write(2,'(A)') 'reset_timestep 0'
    write(2,'(A)') 'dump mydump all custom 1 SREV-unfolded.xyz xu yu zu'
    write(2,'(A)') 'dump_modify mydump sort id header no'
    write(2,'(A)') 'run 0'
    write(2,'(A,F10.3,A)') 'delete_atoms overlap ',2*node_radius,' all all'
    write(2,'(A)') 'shell sleep 0.1'
    write(2,'(A)') 'label finish'
    write(2,'(A)') 'shell rm -f gr-init.txt gr.txt _stop _analyze_gr.pl'
    close(2)

    if (.not.silent) print '(A)','created config.data'
    if (.not.silent) print '(A)',launcher
    call system(launcher)
    ! check if successful
    inquire(file='SREV-unfolded.xyz',EXIST=EX)
    lammps_not_converged = .false.
    if (.not.EX) then
       lammps_not_converged = .true.; print '(A30)','lammps did not converge'
    endif
    inquire(file='_mindistance',exist=EX)
    if (EX) then
       open(2,file='_mindistance'); rewind(2); read(2,*) mindistance; close(2)
       if (mindistance.lt.0.8) then
          lammps_not_converged = .true.; print '(A30)','lammps did not converge. mindistance too small'
       endif
    endif
    call system('rm -f _lammps.in config.data _lammps.in _mindistance')
    if (lammps_not_converged) then
       call system('rm -f SREV-unfolded.xyz')
       error stop 'LAMMPS did not produce a valid SREV configuration'
    else
       if (.not.silent) print '(A)','created SREV-unfolded.xyz'
    endif
    open(2,file='SREV-unfolded.xyz'); rewind(2)
    do id=0,N
       read(2,*) xu(id,:)
    enddo
    close(2)
    ! report mean SRRW-SREV displacement
    mean_SRRW_SREV_displacement = 0.D0
    do id=0,N
       dx = xu(id,:)-xu_init(id,:)
       mean_SRRW_SREV_displacement = mean_SRRW_SREV_displacement + dsqrt(dot_product(dx,dx))
    enddo
    mean_SRRW_SREV_displacement = mean_SRRW_SREV_displacement/dble(N+1)
    if (.not.silent) print '(A30,F10.3)','SRRW-SREV displacement',mean_SRRW_SREV_displacement
    return
  end subroutine create_SREV_via_lammps

  ! ================================================================== 
  function L_for_sphere()
    implicit none
    real*8   :: L_for_sphere
    L_for_sphere = (8.D0*dble(N)*node_radius**3/phi)**(1.D0/3.D0)
    if (.not.silent) print '(A30,F10.3)','alpha =',alpha
    if (.not.silent) print '(A30,F10.3)','phi = ',phi
    if (.not.silent) print '(A30,F10.3)','node radius = ',node_radius
    if (.not.silent) print '(A30,F10.3)','L (sphere) = ',L_for_sphere
    return
  end function L_for_sphere

  ! ================================================================== 
  function L_for_box()
    implicit none
    real*8   :: L_for_box
    real*8, parameter    :: pi = dacos(-1.D0)
    L_for_box = (dble(N)*(4.D0*pi/3.D0)*node_radius**3/phi)**(1.D0/3.D0)
    if (.not.silent) print '(A30,F10.3)','alpha =',alpha
    if (.not.silent) print '(A30,F10.3)','phi = ',phi
    if (.not.silent) print '(A30,F10.3)','node radius = ',node_radius
    if (.not.silent) print '(A30,F10.3)','L (box) = ',L_for_box
    return
  end function L_for_box

  ! =================================================================== 
  subroutine setup_conversion_and_mean_bondlen_for_this_SRRW_model
    real*8   :: sum
    real     :: myrand
    integer  :: trials = 2000000
    integer  :: j
    logical  :: EX
    integer  :: tmp_N
    real*8   :: tmp_alpha,tmp_bmax,tmp_phi
    chi = 1.D0-bmax**(-1.D0-alpha)
    sum = 0.D0
    do j=1,trials
       call random_number(myrand)
       sum = sum + dexp(-dlog(1.D0-chi*myrand)/(1.D0+alpha))
    enddo
    mean_bondlen_for_this_SRRW_model = sum/dble(trials)

    if (first_config) then  ! added 30 jan 2025
       conversion_bondlen = mean_bondlen_for_this_SRRW_model
    else
       inquire(file='_SREV-bondlength.txt',exist=EX)
       if (.not.EX) stop 'you have to set first_config = .true. to create _SREV-bondlength.txt, afterwards use first_config = .false.'
       if (.not.silent) print '(A)','## reading _SREV-bondlength.txt ..'
       open(22,file='_SREV-bondlength.txt'); rewind(22)
       read(22,*) conversion_bondlen
       read(22,*) tmp_N,tmp_alpha,tmp_bmax,tmp_phi
       if (N.ne.tmp_n.or.max(dabs(tmp_alpha-alpha),max(dabs(tmp_bmax-bmax),dabs(tmp_phi-phi))).gt.1.D-3) then
          stop 'you have to set first_config = .true. to create _SREV-bondlength.txt, afterwards use first_config = .false.'
       endif
    endif

    return
  end subroutine setup_conversion_and_mean_bondlen_for_this_SRRW_model

  ! ===================================================================== 
  subroutine get_clock(tictoc)
    real*8  :: tictoc
    integer :: count, count_max, count_rate
    call system_clock(count, count_rate, count_max)
    tictoc = dble(count)/count_rate
    return
  END subroutine get_clock

  ! ===================================================================== generic single-type output
  subroutine save_SREV_data_file
    implicit none
    integer  :: id
    if (.not.silent) print '(A)','subroutine save_SREV_data_file'
    open(2,file='SREV.data'); rewind(2)
    write(2,'(A)') 'SREV 2025'
    write(2,'(A)') ' '
    write(2,'(I10,A)') N+1,' atoms'
    write(2,'(I10,A)') 1,' atom types'
    write(2,'(I10,A)') N,' bonds'
    write(2,'(A)') '1 bond types'
    write(2,'(A)') ' '
    write(2,'(2(F20.8,1x),A)') -L/2.D0,L/2.D0,' xlo xhi'
    write(2,'(2(F20.8,1x),A)') -L/2.D0,L/2.D0,' ylo yhi'
    write(2,'(2(F20.8,1x),A)') -L/2.D0,L/2.D0,' zlo zhi'
    write(2,'(A)') ' '
    write(2,'(A)') 'Atoms # angle'
    write(2,'(A)') ' '
    do id=0,N
       write(2,'(I10,1x,I1,1x,I6,1x,3(F15.8,1x),A)') 1+id,1,1,xu(id,:),' 0 0 0'
    enddo
    write(2,'(A)') ' '
    write(2,'(A)') 'Bonds'
    write(2,'(A)') ' '
    do id=1,N
       write(2,'(I10,1x,I1,1x,2(I6,1x))') id,1,id,id+1
    enddo
    close(2)
    return
  end subroutine save_SREV_data_file

  subroutine HEADER 
    if (silent) return
    print '(A)',' '
    print '(A)','    __   __   __' 
    print '(A)','   |__  |__| |__ \  /   (c) 2026 martin kroger'
    print '(A)','    __| |  \ |__  \/             marcelo carignano'
    print '(A)',' '
  return
  end

end module shared

! ===================================================================== MAIN PROGRAM

program code
  use shared
  implicit none
  real*8   :: tic,toc
  call get_clock(tic)

  call setup_seed_and_constants
  call HEADER
  call setup_conversion_and_mean_bondlen_for_this_SRRW_model    ! creates bp_for_id() and id_for_bp() functions

  ! ---- create SRRW in periodic box
  if (sphere.eq.0) then
    L = L_for_box(); call SRRW_in_periodic_box; if (.not.silent) print '(A)','SREV in box'
  endif

  ! ---- alternative spherical confinement
  if (sphere.eq.1) then
    L = L_for_sphere(); call SRRW_in_sphere; if (.not.silent) print '(A)','SREV in sphere'
  endif

  ! ---- dump the SRRW
  ! call dump_folded_walk(N,xu,type,'SRRW-folded.dump')
  call dump_unfolded_walk(N,xu,type,'SRRW-unfolded.dump')

  ! ---- create the SREV
  call create_SREV_via_lammps

  ! ---- dump the SREV
  call dump_unfolded_walk(N,xu,type,'SREV-unfolded.dump')

  ! ---- save _mean_bond_length and _SREV_average_linker_bp.txt (needed if first_config = .false. is used!)
  call save_mean_bond_length

  ! ---- save _Cp_id.txt 
  call save_Cp_id_fast                  ! creates fort.42 (use this one if create_contact_map = true)
  ! call save_Cp_id_logarithmic         ! creates fort.42 (equi-logarithmic)
  ! The brute-force routine is O(N**2) and is not used in production runs.
  ! call create_contacts_brute_force

  ! --- save SREV as lammps data file    ! added 23 jan 2025
  call save_SREV_data_file               ! choose version here

1 continue

  call get_clock(toc)
  if (.not.silent) print '(A,F6.2,A)','walltime ',toc-tic,' secs'

  stop
end program code

