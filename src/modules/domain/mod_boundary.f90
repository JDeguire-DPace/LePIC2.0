module mod_boundary
  use iso_fortran_env  , only: real64, int32
  use mod_geometry     , only: Domain
  use mod_InputFiles
  use mod_functionsText
  use mod_io_hdf5      , only: save_array_h5, save_array_h5_wo_ghost 
  implicit none
  private

  ! --------------------------------------------------------------------
  ! Boundary conventions (matches legacy generate_boundary):
  !
  !   bcnd(i,j,k):
  !      -2 : Neumann bnd on left side (x=0 ghost cells)
  !      -1 : plasma interior
  !       0 : “inactive” / vacuum / periodic ghosts / hole edges
  !      >0 : wall label (metal or dielectric)
  !
  !   wall_type(label) (final value, always >= 0):
  !      1 : metal
  !      2 : dielectric surface in XZ plane (Y = const)
  !      3 : dielectric surface in YZ plane (first plane, X = const)
  !      4 : dielectric surface in YZ plane (second plane, X = const)
  !
  ! During construction, we use a signed work array wall_type_work:
  !   wall_type_work(label) < 0  => Z-periodic counterpart (legacy)
  ! --------------------------------------------------------------------
  type :: Boundary
    ! Main arrays
    integer(int32), allocatable :: bcnd(:,:,:)                 !! bnd condition indices
    real(real64)  , allocatable :: phi(:,:,:)                  !! electrostatic potential
    integer(int32), allocatable :: wall_type(:)                !! final wall type per label
    real(real64)  , allocatable :: potential_from_boundary(:)  !! label potentials
    integer(int32)              :: nbPotentialLabel = 0        !! number of labels (ngrid)

    ! Legacy-like metadata and flags
    logical                     :: has_cylinders_along_x   = .false.
    logical                     :: has_cylinders_along_z   = .false.
    logical                     :: has_periodic_boundary   = .false.
    logical                     :: has_periodic_boundary_z = .false.
    logical                     :: has_dielectrics         = .false.
    logical                     :: has_neumann_left        = .false.
    logical                     :: has_cylindrical_heating = .false.
    integer(int32)              :: dielectric_plane_index(4) = 0
    integer(int32)              :: num_cathode_surfaces      = 0
    real(real64)                :: secondary_z_positions(2)  = 0.0_real64  !! m

    ! Grid aperture metadata
    real(real64)                :: grid_aperture_span_y = 0.0_real64  !! m
    real(real64)                :: grid_aperture_span_z = 0.0_real64  !! m
    real(real64)                :: grid_aperture_start_x = 0.0_real64 !! m
    integer(int32)              :: grid_aperture_center_index_x = 0
    real(real64)                :: grid_effective_area = 0.0_real64   !! m^2

    ! For diagnostics
    real(real64)                :: last_cylinder_radius_cm = 0.0_real64
  end type Boundary

  public :: Boundary, build_boundary

  contains
  !======================================================================
  ! Public entry point: build_boundary
  !======================================================================
  subroutine build_boundary(bnd, dom, &
                            xl_pow, xr_pow, &
                            use_grid_apertures, grid_label, &
                            aperture_span_y_cm, aperture_span_z_cm, &
                            num_apertures_y, num_apertures_z, &
                            gamma_secondary, injection_option, secondary_grid_label, &
                            mpi_rank)
    type(Boundary), intent(inout) :: bnd
    type(Domain)  , intent(in)    :: dom

    ! Optional inputs (legacy-style)
    real(real64),  intent(in), optional :: xl_pow, xr_pow
    logical,       intent(in), optional :: use_grid_apertures
    integer(int32),intent(in), optional :: grid_label, num_apertures_y, num_apertures_z
    real(real64),  intent(in), optional :: aperture_span_y_cm, aperture_span_z_cm
    real(real64),  intent(in), optional :: gamma_secondary
    integer(int32),intent(in), optional :: injection_option, secondary_grid_label
    integer,       intent(in), optional :: mpi_rank

    ! Domain grid size
    integer :: num_cells_x, num_cells_y, num_cells_z

    ! Geometry (cm)
    real(real64) :: box_length_x_cm, box_length_y_cm, box_length_z_cm
    real(real64) :: cell_size_cm(3)   ! [dx,dy,dz] in cm

    ! Options (local copies)
    logical        :: enable_grid_apertures
    integer(int32) :: grid_label_local, num_ap_y_local, num_ap_z_local
    integer(int32) :: injection_option_local, secondary_grid_label_local
    real(real64)   :: aperture_span_y_local_cm, aperture_span_z_local_cm
    real(real64)   :: gamma_secondary_local
    real(real64)   :: xl_pow_local, xr_pow_local
    integer        :: mpi_rank_local

    ! Flags (integer form, like legacy)
    integer :: cylinder_flag_from_header_x  ! from sign of ymax
    integer :: cylinder_flag_z              ! from sign of zmax
    integer :: periodic_flag_any
    integer :: periodic_flag_z
    integer :: neumann_left_flag
    integer :: dielectric_present_flag
    integer :: yz_dielectric_plane_count

    ! I/O
    character(len=256)  :: geometry_file, boundary_file
    character(len=1024) :: text_line
    integer :: unit_geom, ios

    ! Single wall row
    real(real64) :: wall_x_min_cm, wall_y_min_cm, wall_z_min_cm
    real(real64) :: wall_x_max_cm, wall_y_max_cm, wall_z_max_cm
    integer      :: wall_label

    ! Work array for wall types
    integer(int32), allocatable :: wall_type_work(:)

    ! Grid segment indices for apertures
    integer :: grid_index_x_min, grid_index_x_max

    ! Misc
    integer      :: label_index
    real(real64) :: grid_area_cm2, pi

    pi = acos(-1.0_real64)

    ! --- dom sizes ---
    num_cells_x = dom%ncell_x1
    num_cells_y = dom%ncell_x2
    num_cells_z = dom%ncell_x3

    bnd%nbPotentialLabel = dom%get_nbPotentialLabel()
    if (bnd%nbPotentialLabel <= 0) error stop 'build_boundary: invalid nbPotentialLabel'

    ! --- optionals -> locals ---
    enable_grid_apertures       = .false.; if (present(use_grid_apertures))       enable_grid_apertures       = use_grid_apertures
    grid_label_local            = 0      ; if (present(grid_label))              grid_label_local            = grid_label
    aperture_span_y_local_cm    = 0.0d0  ; if (present(aperture_span_y_cm))      aperture_span_y_local_cm    = aperture_span_y_cm
    aperture_span_z_local_cm    = 0.0d0  ; if (present(aperture_span_z_cm))      aperture_span_z_local_cm    = aperture_span_z_cm
    num_ap_y_local              = 1      ; if (present(num_apertures_y))         num_ap_y_local              = max(1, num_apertures_y)
    num_ap_z_local              = 1      ; if (present(num_apertures_z))         num_ap_z_local              = max(1, num_apertures_z)
    gamma_secondary_local       = 0.0d0  ; if (present(gamma_secondary))         gamma_secondary_local       = gamma_secondary
    injection_option_local      = 0      ; if (present(injection_option))        injection_option_local      = injection_option
    secondary_grid_label_local  = 0      ; if (present(secondary_grid_label))    secondary_grid_label_local  = secondary_grid_label
    xl_pow_local                = -1.0d300; if (present(xl_pow))                 xl_pow_local                = xl_pow
    xr_pow_local                =  1.0d300; if (present(xr_pow))                 xr_pow_local                = xr_pow
    mpi_rank_local              = 0      ; if (present(mpi_rank))                mpi_rank_local              = mpi_rank

    ! --- init arrays + flags ---
    call initialize_boundary_arrays_and_flags(bnd, num_cells_x, num_cells_y, num_cells_z)

    ! --- read V(label) from bnd.inp ---
    boundary_file = trim(find_boundary_file())
    call read_label_potentials(boundary_file, bnd%nbPotentialLabel, bnd%potential_from_boundary)

    ! --- read geometry header ---
    geometry_file = trim(find_geometry_file())
    open(newunit=unit_geom, file=geometry_file, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'build_boundary: cannot open '//geometry_file

    ! skip first line: ncells
    call next_data_line(unit_geom, ios); read(unit_geom,*)

    ! second line: box size (cm with signs)
    call next_data_line(unit_geom, ios)
    read(unit_geom,*,iostat=ios) box_length_x_cm, box_length_y_cm, box_length_z_cm
    if (ios /= 0) error stop 'build_boundary: bad box-size line in geometry.inp'

    ! Neumann left (x<0)
    neumann_left_flag = 0
    if (box_length_x_cm < 0.0_real64) then
      neumann_left_flag = 1
      box_length_x_cm   = abs(box_length_x_cm)
    end if

    ! cylindrical along x from sign(ymax)
    cylinder_flag_from_header_x = 0
    if (box_length_y_cm < 0.0_real64) cylinder_flag_from_header_x = 1
    box_length_y_cm = abs(box_length_y_cm)

    ! cylindrical along z from sign(zmax)
    cylinder_flag_z = 0
    if (box_length_z_cm < 0.0_real64) cylinder_flag_z = 1
    box_length_z_cm = abs(box_length_z_cm)

    cell_size_cm(1) = box_length_x_cm / real(num_cells_x, real64)
    cell_size_cm(2) = box_length_y_cm / real(num_cells_y, real64)
    cell_size_cm(3) = box_length_z_cm / real(num_cells_z, real64)

    periodic_flag_any        = 0
    periodic_flag_z          = 0
    dielectric_present_flag  = 0
    yz_dielectric_plane_count= 0
    bnd%dielectric_plane_index = 0
    bnd%secondary_z_positions  = 0.0_real64
    bnd%num_cathode_surfaces   = 0
    bnd%grid_aperture_span_y   = 0.0_real64
    bnd%grid_aperture_span_z   = 0.0_real64
    bnd%grid_aperture_start_x  = 0.0_real64
    bnd%grid_aperture_center_index_x = 0
    bnd%grid_effective_area    = 0.0_real64
    bnd%last_cylinder_radius_cm= 0.0_real64
    grid_index_x_min = 0
    grid_index_x_max = 0

    allocate(wall_type_work(0:bnd%nbPotentialLabel))
    wall_type_work = 0

    ! --- main loop over wall segments ---
    do
      call next_data_line(unit_geom, ios)
      if (ios /= 0) exit

      read(unit_geom,'(A)',iostat=ios) text_line
      if (ios /= 0) exit
      call strip_inline_comment(text_line)
      if (len_trim(text_line) == 0) cycle

      if (is_single_integer(text_line)) exit  ! trailer

      ios = 0
      read(text_line,*,iostat=ios) &
           wall_x_min_cm, wall_y_min_cm, wall_z_min_cm, &
           wall_x_max_cm, wall_y_max_cm, wall_z_max_cm, wall_label
      if (ios /= 0) error stop 'build_boundary: malformed wall row in geometry.inp'

      call process_wall_segment( &
           wall_x_min_cm, wall_y_min_cm, wall_z_min_cm, &
           wall_x_max_cm, wall_y_max_cm, wall_z_max_cm, wall_label, &
           num_cells_x, num_cells_y, num_cells_z, cell_size_cm, &
           box_length_x_cm, box_length_y_cm, box_length_z_cm, &
           cylinder_flag_from_header_x, cylinder_flag_z, &
           periodic_flag_any, periodic_flag_z, dielectric_present_flag, &
           yz_dielectric_plane_count, &
           bnd, wall_type_work, &
           xl_pow_local, xr_pow_local, &
           grid_label_local, gamma_secondary_local, injection_option_local, &
           secondary_grid_label_local, mpi_rank_local, &
           grid_index_x_min, grid_index_x_max )
    end do

    close(unit_geom)

    ! --- Neumann at left ---
    call apply_neumann_left(bnd, neumann_left_flag)

    ! --- periodic Z caps ---
    call apply_periodic_z_caps(bnd, wall_type_work, periodic_flag_any, &
                               num_cells_x, num_cells_y, num_cells_z)

    ! finalize wall_type = |work|
    do label_index = 1, bnd%nbPotentialLabel
      bnd%wall_type(label_index) = abs(wall_type_work(label_index))
    end do

    ! --- grid apertures ---
    grid_area_cm2 = 0.0_real64
    if (enable_grid_apertures) then
      call draw_grid_apertures(bnd, cell_size_cm, &
                               box_length_x_cm, box_length_y_cm, box_length_z_cm, &
                               aperture_span_y_local_cm, aperture_span_z_local_cm, &
                               num_ap_y_local, num_ap_z_local, grid_label_local, &
                               mpi_rank_local, grid_index_x_min, grid_index_x_max, &
                               grid_area_cm2)
      bnd%grid_effective_area = grid_area_cm2 * 1.0d-4
    end if

    ! --- convert cm → m, set logical flags ---
    call convert_units_and_set_logical_flags(bnd, cell_size_cm, &
                                             box_length_x_cm, box_length_y_cm, box_length_z_cm, &
                                             cylinder_flag_from_header_x, cylinder_flag_z, &
                                             periodic_flag_any, periodic_flag_z, &
                                             dielectric_present_flag, neumann_left_flag)

    ! --- save to HDF5 ---
    call write_boundary_hdf5(bnd, dom)

    deallocate(wall_type_work)
  end subroutine build_boundary



  !======================================================================
  subroutine initialize_boundary_arrays_and_flags(bnd, nx, ny, nz)
    type(Boundary), intent(inout) :: bnd
    integer,        intent(in)    :: nx, ny, nz

    if (allocated(bnd%bcnd)) deallocate(bnd%bcnd, bnd%phi)
    allocate(bnd%bcnd(0:nx+2, 0:ny+2, 0:nz+2))
    allocate(bnd%phi (0:nx+2, 0:ny+2, 0:nz+2))

    bnd%bcnd = 0_int32
    bnd%phi  = 0.0_real64

    if (allocated(bnd%wall_type)) deallocate(bnd%wall_type)
    allocate(bnd%wall_type(0:bnd%nbPotentialLabel))
    bnd%wall_type = 0_int32

    if (allocated(bnd%potential_from_boundary)) deallocate(bnd%potential_from_boundary)
    allocate(bnd%potential_from_boundary(bnd%nbPotentialLabel))

    bnd%has_cylinders_along_x   = .false.
    bnd%has_cylinders_along_z   = .false.
    bnd%has_periodic_boundary   = .false.
    bnd%has_periodic_boundary_z = .false.
    bnd%has_dielectrics         = .false.
    bnd%has_neumann_left        = .false.
    bnd%has_cylindrical_heating = .false.
    bnd%dielectric_plane_index  = 0
    bnd%num_cathode_surfaces    = 0
    bnd%secondary_z_positions   = 0.0_real64
    bnd%grid_aperture_span_y    = 0.0_real64
    bnd%grid_aperture_span_z    = 0.0_real64
    bnd%grid_aperture_start_x   = 0.0_real64
    bnd%grid_aperture_center_index_x = 0
    bnd%grid_effective_area     = 0.0_real64
    bnd%last_cylinder_radius_cm = 0.0_real64
  end subroutine initialize_boundary_arrays_and_flags



  !======================================================================
  subroutine read_label_potentials(boundary_file, num_labels, potentials)
    character(len=*), intent(in)    :: boundary_file
    integer,          intent(in)    :: num_labels
    real(real64),     intent(inout) :: potentials(num_labels)

    integer :: unit_id, ios, i

    open(newunit=unit_id, file=boundary_file, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'build_boundary: cannot open '//boundary_file

    do i = 1, num_labels
      call next_data_line(unit_id, ios)
      if (ios /= 0) error stop 'build_boundary: missing potentials in bnd.inp'
      read(unit_id,*,iostat=ios) potentials(i)
      if (ios /= 0) error stop 'build_boundary: bad potential line in bnd.inp'
    end do

    close(unit_id)
  end subroutine read_label_potentials



  !======================================================================
  subroutine process_wall_segment(x_left_cm,  y_left_cm,  z_left_cm, &
                                  x_right_cm, y_right_cm, z_right_cm, wall_label, &
                                  nx, ny, nz, cell_size_cm, &
                                  box_x_cm, box_y_cm, box_z_cm, &
                                  cylinder_flag_header_x, cylinder_flag_z, &
                                  periodic_flag_any, periodic_flag_z, dielectric_flag, &
                                  yz_dielectric_count, &
                                  bnd, wall_type_work, &
                                  xl_pow, xr_pow, &
                                  grid_label, gamma_secondary, injection_option, secondary_grid_label, &
                                  mpi_rank, &
                                  grid_ix_min, grid_ix_max)
    real(real64), intent(inout) :: x_left_cm,  y_left_cm,  z_left_cm
    real(real64), intent(inout) :: x_right_cm, y_right_cm, z_right_cm
    integer,      intent(in)    :: wall_label
    integer,      intent(in)    :: nx, ny, nz
    real(real64), intent(in)    :: cell_size_cm(3)
    real(real64), intent(in)    :: box_x_cm, box_y_cm, box_z_cm
    integer,      intent(in)    :: cylinder_flag_header_x, cylinder_flag_z
    integer,      intent(inout) :: periodic_flag_any, periodic_flag_z, dielectric_flag
    integer,      intent(inout) :: yz_dielectric_count
    type(Boundary), intent(inout) :: bnd
    integer(int32), intent(inout) :: wall_type_work(0:bnd%nbPotentialLabel)
    real(real64),   intent(in)    :: xl_pow, xr_pow, gamma_secondary
    integer,        intent(in)    :: grid_label, injection_option, secondary_grid_label, mpi_rank
    integer,        intent(inout) :: grid_ix_min, grid_ix_max

    logical      :: use_cylinder_x, use_cylinder_z
    real(real64) :: R_cm
    integer      :: label_abs
    integer      :: ixL, ixR, iyL, iyR, izL, izR
    real(real64) :: xL0, yL0, zL0, xR0, yR0, zR0

    xL0 = x_left_cm ; yL0 = y_left_cm ; zL0 = z_left_cm
    xR0 = x_right_cm; yR0 = y_right_cm; zR0 = z_right_cm

    label_abs = abs(wall_label)

    ! default to metal
    if (label_abs > 0 .and. label_abs <= bnd%nbPotentialLabel) then
      if (wall_type_work(label_abs) == 0) wall_type_work(label_abs) = 1
    end if

    ! --- cylinder flags for this wall ---
    use_cylinder_x = (cylinder_flag_header_x == 1)
    use_cylinder_z = (cylinder_flag_z        == 1)

    if (y_left_cm < 0.0_real64) then
      if (cylinder_flag_header_x == 1) then
        use_cylinder_x = .false.
      else
        use_cylinder_x = .true.
      end if
    end if
    y_left_cm = abs(y_left_cm)

    ! radius
    R_cm = 0.0_real64
    if (use_cylinder_x) then
      R_cm = min( (y_right_cm - y_left_cm)/2.0d0, (z_right_cm - z_left_cm)/2.0d0 )
      if (x_left_cm <= xl_pow .and. x_right_cm >= xr_pow) bnd%has_cylindrical_heating = .true.
    end if
    if (use_cylinder_z) then
      R_cm = min( (x_right_cm - x_left_cm)/2.0d0, (y_right_cm - y_left_cm)/2.0d0 )
    end if
    bnd%last_cylinder_radius_cm = R_cm

    ! periodic label
    if (wall_label == 0) periodic_flag_any = 1

    ! dielectric YZ planes (X=const)
    if (x_right_cm < 0.0_real64) then
      x_right_cm = abs(x_right_cm)
      yz_dielectric_count = yz_dielectric_count + 1
      if (yz_dielectric_count == 1) then
        wall_type_work(label_abs) = 3
      else
        wall_type_work(label_abs) = 4
      end if
      dielectric_flag = 1
      bnd%dielectric_plane_index(wall_type_work(label_abs)) = int( x_left_cm / cell_size_cm(1) ) + 1
    end if

    ! dielectric XZ planes (Y=const)
    if (y_right_cm < 0.0_real64) then
      y_right_cm = abs(y_right_cm)
      wall_type_work(label_abs) = 2
      dielectric_flag           = 1
      bnd%dielectric_plane_index(2) = int( y_left_cm / cell_size_cm(2) ) + 1
    end if

    ! Z-periodic
    if (z_right_cm < 0.0_real64) then
      z_right_cm = abs(z_right_cm)
      wall_type_work(label_abs) = -wall_type_work(label_abs)
      periodic_flag_any = 1
      periodic_flag_z   = 1
    end if

    ! label safety
    if (label_abs > bnd%nbPotentialLabel) then
      if (mpi_rank == 0) then
        write(*,*) 'Insufficient wall labels; correct bnd.inp'
      end if
      error stop 'build_boundary: wall label exceeds nbPotentialLabel'
    end if

    ! indices 0..n+2
    ixL = int( x_left_cm  / cell_size_cm(1) ) + 1; if (ixL <= 1) ixL = 0
    ixR = int( x_right_cm / cell_size_cm(1) ) + 1; if (ixR >= nx) ixR = nx+2
    iyL = int( y_left_cm  / cell_size_cm(2) ) + 1; if (iyL <= 1) iyL = 0
    iyR = int( y_right_cm / cell_size_cm(2) ) + 1; if (iyR >= ny) iyR = ny+2
    izL = int( z_left_cm  / cell_size_cm(3) ) + 1; if (izL <= 1) izL = 0
    izR = int( z_right_cm / cell_size_cm(3) ) + 1; if (izR >= nz) izR = nz+2

    ! secondary emission meta
    if ( (gamma_secondary > 0.0d0 .or. abs(injection_option) == 4) .and. secondary_grid_label == label_abs ) then
      if (zR0 /= zL0) then
        write(*,*) 'Warning: secondary emission only along Oz; correct geometry.'
        error stop 'build_boundary: invalid secondary emission plane'
      end if
      if (zL0 == 0.0d0) then
        bnd%secondary_z_positions(1) = zR0
        bnd%num_cathode_surfaces     = bnd%num_cathode_surfaces + 1
      end if
      if (zR0 == box_z_cm) then
        bnd%secondary_z_positions(2) = zL0
        bnd%num_cathode_surfaces     = bnd%num_cathode_surfaces + 1
      end if
    end if

    ! grid region
    if (label_abs == grid_label) then
      grid_ix_min  = ixL
      grid_ix_max  = ixR
      bnd%grid_aperture_span_y = min( (iyR - iyL)*cell_size_cm(2), box_y_cm )
      bnd%grid_aperture_span_z = min( (izR - izL)*cell_size_cm(3), box_z_cm )
      bnd%grid_aperture_center_index_x = nint(real(grid_ix_min + grid_ix_max, real64)/2.0d0)
      bnd%grid_aperture_start_x        = (ixL-1)*cell_size_cm(1)
    end if

    ! periodic region in X
    if (wall_label == 0) then
      ixL = max(ixL, 2)
      ixR = min(ixR, nx)
    end if

    ! first pass: fill wall block
    if (wall_label >= 0) call fill_wall_block(bnd, ixL, ixR, wall_label)

    ! restrict to interior for carving
    if (wall_label >= 0) then
      ixL = max(ixL, 2); ixR = min(ixR, nx)
      iyL = max(iyL, 2); iyR = min(iyR, ny)
      izL = max(izL, 2); izR = min(izR, nz)
    end if

    call carve_interior_with_cylinders(bnd, wall_label, &
                                       ixL, ixR, iyL, iyR, izL, izR, &
                                       cell_size_cm, box_x_cm, box_y_cm, box_z_cm, &
                                       R_cm, use_cylinder_x, (cylinder_flag_z == 1))
  end subroutine process_wall_segment



  !======================================================================
  subroutine fill_wall_block(bnd, ixL, ixR, wall_label)
    type(Boundary), intent(inout) :: bnd
    integer,        intent(in)    :: ixL, ixR, wall_label

    integer :: ix, iy, iz
    integer :: jmin, jmax, kmin, kmax

    if (ixL > ixR) return

    jmin = lbound(bnd%bcnd,2); jmax = ubound(bnd%bcnd,2)
    kmin = lbound(bnd%bcnd,3); kmax = ubound(bnd%bcnd,3)

    do iz = kmin, kmax
      do iy = jmin, jmax
        do ix = ixL, ixR
          bnd%bcnd(ix,iy,iz) = wall_label
          if (wall_label > 0) bnd%phi(ix,iy,iz) = bnd%potential_from_boundary(wall_label)
        end do
      end do
    end do
  end subroutine fill_wall_block



  !======================================================================
  subroutine carve_interior_with_cylinders(bnd, wall_label, &
                                           ixL, ixR, iyL, iyR, izL, izR, &
                                           cell_size_cm, box_x_cm, box_y_cm, box_z_cm, &
                                           radius_cm, use_cyl_x, use_cyl_z)
    type(Boundary), intent(inout) :: bnd
    integer,        intent(in)    :: wall_label
    integer,        intent(in)    :: ixL, ixR, iyL, iyR, izL, izR
    real(real64),   intent(in)    :: cell_size_cm(3)
    real(real64),   intent(in)    :: box_x_cm, box_y_cm, box_z_cm
    real(real64),   intent(in)    :: radius_cm
    logical,        intent(in)    :: use_cyl_x, use_cyl_z

    integer :: ix, iy, iz
    real(real64) :: x_cm, y_cm, z_cm

    if (ixL > ixR .or. iyL > iyR .or. izL > izR) return

    do iz = izL, izR
      z_cm = (iz-1)*cell_size_cm(3)
      do iy = iyL, iyR
        y_cm = (iy-1)*cell_size_cm(2)
        do ix = ixL, ixR
          x_cm = (ix-1)*cell_size_cm(1)

          if (wall_label >= 0) then
            if (use_cyl_z .and. wall_label > 0) then
              if (( (x_cm - box_x_cm/2.0d0)**2 + (y_cm - box_y_cm/2.0d0)**2 ) > radius_cm**2) cycle
            end if
            if (use_cyl_x .and. wall_label > 0) then
              if (( (y_cm - box_y_cm/2.0d0)**2 + (z_cm - box_z_cm/2.0d0)**2 ) > radius_cm**2) cycle
            end if

            bnd%bcnd(ix,iy,iz) = -1
            bnd%phi (ix,iy,iz) = 0.0d0

          else
            if (use_cyl_z) then
              if (( (x_cm - box_x_cm/2.0d0)**2 + (y_cm - box_y_cm/2.0d0)**2 ) > radius_cm**2) cycle
            end if
            if (use_cyl_x) then
              if (( (y_cm - box_y_cm/2.0d0)**2 + (z_cm - box_z_cm/2.0d0)**2 ) > radius_cm**2) cycle
            end if

            bnd%bcnd(ix,iy,iz) = abs(wall_label)
            if (abs(wall_label) > 0) bnd%phi(ix,iy,iz) = bnd%potential_from_boundary(abs(wall_label))
          end if
        end do
      end do
    end do
  end subroutine carve_interior_with_cylinders



  !======================================================================
  subroutine apply_neumann_left(bnd, neumann_flag)
    type(Boundary), intent(inout) :: bnd
    integer,        intent(in)    :: neumann_flag

    if (neumann_flag /= 1) return
    bnd%bcnd(0:1,:,:) = -2
    bnd%phi (0:1,:,:) = 0.0d0
  end subroutine apply_neumann_left



  !======================================================================
  subroutine apply_periodic_z_caps(bnd, wall_type_work, periodic_flag_any, &
                                   nx, ny, nz)
    type(Boundary), intent(inout) :: bnd
    integer(int32), intent(in)    :: wall_type_work(:)
    integer,        intent(in)    :: periodic_flag_any, nx, ny, nz

    integer :: ix, iy, label_here

    if (periodic_flag_any == 0) return

    do iy = 0, ny+2
      do ix = 0, nx+2
        label_here = bnd%bcnd(ix,iy,1)
        if (label_here < 0) cycle
        if (label_here > ubound(wall_type_work,1)) cycle
        if (wall_type_work(label_here) < 0) then
          if (bnd%bcnd(ix,iy,2) == -1) then
            bnd%bcnd(ix,iy,0:1)               = 0
            bnd%bcnd(ix,iy,nz+1:nz+2)         = 0
            bnd%phi (ix,iy,0:1)               = 0.0d0
            bnd%phi (ix,iy,nz+1:nz+2)         = 0.0d0
          end if
        end if
      end do
    end do
  end subroutine apply_periodic_z_caps



  !======================================================================
  subroutine draw_grid_apertures(bnd, cell_size_cm, &
                                 box_x_cm, box_y_cm, box_z_cm, &
                                 Lhy_cm, Lhz_cm, &
                                 ny_ap, nz_ap, grid_label, &
                                 mpi_rank, grid_ix_min, grid_ix_max, &
                                 grid_area_cm2)
    type(Boundary), intent(inout) :: bnd
    real(real64) , intent(in)     :: cell_size_cm(3)
    real(real64) , intent(inout)  :: box_x_cm, box_y_cm, box_z_cm
    real(real64) , intent(inout)  :: Lhy_cm, Lhz_cm
    integer      , intent(in)     :: ny_ap, nz_ap, grid_label
    integer      , intent(in)     :: mpi_rank
    integer      , intent(in)     :: grid_ix_min, grid_ix_max
    real(real64) , intent(out)    :: grid_area_cm2

    integer :: num_cells_x, num_cells_y, num_cells_z
    integer :: half_cells_y, half_cells_z
    integer :: cnt_y, cnt_z
    integer :: iyL, iyR, izL, izR
    integer :: ix, iy, iz, value_inside
    real(real64) :: yd_cm, zd_cm, y_cm, z_cm
    real(real64) :: R_cm, aperture_area_cm2, pi

    pi = acos(-1.0_real64)

    num_cells_x = ubound(bnd%bcnd,1) - 2
    num_cells_y = ubound(bnd%bcnd,2) - 2
    num_cells_z = ubound(bnd%bcnd,3) - 2

    grid_area_cm2 = 0.0_real64
    if (grid_label <= 0 .or. grid_ix_min == 0 .or. grid_ix_max == 0) return

    half_cells_y = nint( Lhy_cm / (2.0d0*cell_size_cm(2)) ) + 1
    half_cells_z = nint( Lhz_cm / (2.0d0*cell_size_cm(3)) ) + 1

    if (Lhy_cm > box_y_cm) Lhy_cm = box_y_cm
    if (Lhz_cm > box_z_cm) Lhz_cm = box_z_cm

    if (Lhy_cm == Lhz_cm) then
      R_cm            = Lhy_cm / 2.0d0
      aperture_area_cm2 = pi*R_cm*R_cm
    else
      R_cm            = 0.0d0
      aperture_area_cm2 = Lhy_cm*Lhz_cm
    end if

    if (mpi_rank == 0) then
      write(*,'(1x,"Corrected dimensions of holes, Ly(cm)= ",f6.2,", Lz(cm)= ",f6.2)') &
           (2*half_cells_y+1)*cell_size_cm(2), (2*half_cells_z+1)*cell_size_cm(3)
    end if

    do cnt_z = 1, nz_ap
      zd_cm = -(nz_ap-1)*3.0d0*Lhz_cm/4.0d0 + (cnt_z-1)*3.0d0*Lhz_cm/2.0d0 + box_z_cm/2.0d0
      izL = nint( zd_cm / cell_size_cm(3) ) + 1 - half_cells_z
      izR = izL + 2*half_cells_z
      if ( (izL < 1 .or. izR > num_cells_z+1) .and. nz_ap > 1 ) then
        if (mpi_rank == 0) then
          write(*,*) 'Warning: too many apertures along Oz'
        end if
        error stop 'build_boundary: too many apertures along Oz'
      end if
      if (izL < 1) izL = 1
      if (izR > num_cells_z+1) izR = num_cells_z+1

      do cnt_y = 1, ny_ap
        yd_cm = -(ny_ap-1)*3.0d0*Lhy_cm/4.0d0 + (cnt_y-1)*3.0d0*Lhy_cm/2.0d0 + box_y_cm/2.0d0
        iyL = nint( yd_cm / cell_size_cm(2) ) + 1 - half_cells_y
        iyR = iyL + 2*half_cells_y
        if ( (iyL < 1 .or. iyR > num_cells_y+1) .and. ny_ap > 1 ) then
          if (mpi_rank == 0) then
            write(*,*) 'Warning: too many apertures along Oy'
          end if
          error stop 'build_boundary: too many apertures along Oy'
        end if
        if (iyL < 1) iyL = 1
        if (iyR > num_cells_y+1) iyR = num_cells_y+1

        do iz = izL, izR
          z_cm = (iz-1)*cell_size_cm(3)
          do iy = iyL, iyR
            y_cm = (iy-1)*cell_size_cm(2)
            value_inside = -1
            if (iy == 1 .or. iy == num_cells_y+1 .or. iz == 1 .or. iz == num_cells_z+1) value_inside = 0

            do ix = grid_ix_min, grid_ix_max
              if (Lhy_cm == Lhz_cm) then
                if ( ( (y_cm-yd_cm)**2 + (z_cm-zd_cm)**2 ) <= R_cm**2 ) then
                  bnd%bcnd(ix,iy,iz) = value_inside
                  bnd%phi (ix,iy,iz) = 0.0d0
                end if
              else
                bnd%bcnd(ix,iy,iz) = value_inside
                bnd%phi (ix,iy,iz) = 0.0d0
              end if
            end do
          end do
        end do
      end do
    end do

    grid_area_cm2 = Lhy_cm*Lhz_cm - &
                    real((ny_ap-1)*(nz_ap-1),real64)*aperture_area_cm2

    if (mpi_rank == 0) then
      write(*,'(1x,"Type of wall surfaces (1=metal, 2 & 3= dielectric):",10(1x,i2))') &
           bnd%wall_type(1:bnd%nbPotentialLabel)
      write(*,'(1x,"Dimension of grid #",i2," : Lgy(cm)= ",f6.2,", Lgz(cm)= ",f6.2)') &
           grid_label, Lhy_cm, Lhz_cm
      write(*,'(1x,"Grid surface (cm2): ",f6.2,", surface occupied by apertures (cm2): ",f6.2)') &
           Lhy_cm*Lhz_cm, real((ny_ap-1)*(nz_ap-1),real64)*aperture_area_cm2
    end if
  end subroutine draw_grid_apertures



  !======================================================================
  subroutine convert_units_and_set_logical_flags(bnd, cell_size_cm, &
                                                 box_x_cm, box_y_cm, box_z_cm, &
                                                 cylinder_flag_header_x, cylinder_flag_z, &
                                                 periodic_flag_any, periodic_flag_z, &
                                                 dielectric_flag, neumann_flag)
    type(Boundary), intent(inout) :: bnd
    real(real64) , intent(inout)  :: cell_size_cm(3)
    real(real64) , intent(inout)  :: box_x_cm, box_y_cm, box_z_cm
    integer      , intent(in)     :: cylinder_flag_header_x, cylinder_flag_z
    integer      , intent(in)     :: periodic_flag_any, periodic_flag_z
    integer      , intent(in)     :: dielectric_flag, neumann_flag

    cell_size_cm                = cell_size_cm*1.0d-2
    box_x_cm                    = box_x_cm*1.0d-2
    box_y_cm                    = box_y_cm*1.0d-2
    box_z_cm                    = box_z_cm*1.0d-2
    bnd%grid_aperture_start_x = bnd%grid_aperture_start_x*1.0d-2
    bnd%grid_aperture_span_y  = bnd%grid_aperture_span_y*1.0d-2
    bnd%grid_aperture_span_z  = bnd%grid_aperture_span_z*1.0d-2
    bnd%secondary_z_positions = bnd%secondary_z_positions*1.0d-2

    bnd%has_cylinders_along_x   = (cylinder_flag_header_x /= 0)
    bnd%has_cylinders_along_z   = (cylinder_flag_z        /= 0)
    bnd%has_periodic_boundary   = (periodic_flag_any      /= 0)
    bnd%has_periodic_boundary_z = (periodic_flag_z        /= 0)
    bnd%has_dielectrics         = (dielectric_flag        /= 0)
    bnd%has_neumann_left        = (neumann_flag           /= 0)
  end subroutine convert_units_and_set_logical_flags



  !======================================================================
  subroutine write_boundary_hdf5(bnd, dom)
    type(Boundary), intent(in) :: bnd
    type(Domain)  , intent(in) :: dom

    integer(int32) :: nghost
    real(real64)   :: origin(3), spacing(3)

    origin  = [0.0_real64, 0.0_real64, 0.0_real64]
    spacing = [dom%dx1, dom%dx2, dom%dx3]
    nghost  = 2

    call save_array_h5_wo_ghost("./Outputs/state.h5","bcnd", bnd%bcnd, origin, spacing, nghost=1, replace=.true.)
    call save_array_h5_wo_ghost("./Outputs/state.h5","phi" , bnd%phi , origin, spacing, nghost=1, replace=.true.)

    call save_array_h5("./Outputs/state.h5","wall_type",               bnd%wall_type,               replace=.true.)
    call save_array_h5("./Outputs/state.h5","potential_from_boundary", bnd%potential_from_boundary, replace=.true.)
  end subroutine write_boundary_hdf5

end module mod_boundary
