module mod_boundary
  use iso_fortran_env, only: real64, int32
  use mod_geometry   , only: Domain
  use mod_InputFiles
  use mod_io_hdf5    , only: save_array_h5
  use mod_functionsText
  implicit none
  private

  type :: Boundary
    integer(int32), allocatable :: bcnd(:,:,:)   !! >0 metal label, 0 interior, -1 non-metal surface
    real(real64)  , allocatable :: phi(:,:,:)    !! Dirichlet potentials on metal labels
    integer(int32), allocatable :: wall_type(:)      !! 1=metal, 2(Y), 3/4(YZ). <0 => same type but periodic in Z
    real(real64)  , allocatable :: potential_from_boundary(:)          !! potentials potential_from_boundary(1..nbPotentialLabel) from boundary.inp
    integer(int32)              :: nbPotentialLabel = 0
    logical                     :: flag_cylinder = .false.
  end type Boundary

  public :: Boundary, build_boundary

  contains

  subroutine build_boundary(bnd, dom)
    type(Boundary), intent(inout) :: bnd
    type(Domain)  , intent(in)    :: dom

    integer :: iu, ios
    integer :: nb_cell_x1,nb_cell_x2,nb_cell_x3, iteration_label, label_potential
    character(len=256) :: geom_file, bound_file
    character(len=1024) :: line
    real(real64) :: x1_left,x2_left,x3_left,x1_right,x2_right,x3_right
    integer :: index_x1_left,index_x1_right,index_x2_left,index_x2_right,index_x3_left,index_x3_right
    real(real64) :: dx,dy,dz, xmax,ymax,zmax
    integer :: max_label, cnt_yz
    integer(int32) :: nghost
    real(real64) :: origin(3), spacing(3)

    nb_cell_x1 = dom%ncell_x1; nb_cell_x2 = dom%ncell_x2; nb_cell_x3 = dom%ncell_x3
    dx = dom%dx1;      dy = dom%dx2;      dz = dom%dx3
    xmax = dom%Length_x1; ymax = dom%Length_x2; zmax = dom%Length_x3

    ! ---- read nbPotentialLabel (last integer line in geometry) ----
    bnd%nbPotentialLabel = dom%get_nbPotentialLabel()
    if (bnd%nbPotentialLabel <= 0) error stop 'build_boundary: invalid nbPotentialLabel'

    ! ---- read potentials potential_from_boundary(1..nbPotentialLabel) ----
    bound_file = trim(find_boundary_file())
    if (allocated(bnd%potential_from_boundary)) deallocate(bnd%potential_from_boundary)
    allocate(bnd%potential_from_boundary(bnd%nbPotentialLabel))
    open(newunit=iu, file=bound_file, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'build_boundary: cannot open '//bound_file
    do iteration_label = 1, bnd%nbPotentialLabel
      call next_data_line(iu, ios); if (ios /= 0) error stop 'build_boundary: missing potentials'
      read(iu,*,iostat=ios) bnd%potential_from_boundary(iteration_label)
      if (ios /= 0) error stop 'build_boundary: bad potential line in boundary.inp'
    end do
    close(iu)

    ! ---- allocate arrays ----
    if (allocated(bnd%bcnd)) deallocate(bnd%bcnd, bnd%phi)
    allocate(bnd%bcnd(0:nb_cell_x1+2,0:nb_cell_x2+2,0:nb_cell_x3+2)); bnd%bcnd = 0
    allocate(bnd%phi (0:nb_cell_x1+2,0:nb_cell_x2+2,0:nb_cell_x3+2)); bnd%phi  = 0.0_real64

    if (allocated(bnd%wall_type)) deallocate(bnd%wall_type)
    allocate(bnd%wall_type(0:bnd%nbPotentialLabel)); bnd%wall_type = 0  ! wall_type(0) unused

    bnd%flag_cylinder = .false.
    max_label = 0
    cnt_yz    = 0        ! to distinguish first/second YZ plane -> 3 or 4

    ! ---- parse geometry wall rows (skip first 2 data lines) ----
    geom_file = trim(find_geometry_file())
    open(newunit=iu, file=geom_file, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'build_boundary: cannot open '//geom_file

    call next_data_line(iu, ios); read(iu,*)
    call next_data_line(iu, ios); read(iu,*)

    do
      call next_data_line(iu, ios); if (ios /= 0) exit
      read(iu,'(A)',iostat=ios) line
      if (ios /= 0) exit
      call strip_inline_comment(line)
      if (len_trim(line) == 0) cycle

      ! --- sentinel check: a line with a single integer => nbPotentialLabel at end
      if (is_single_integer(line) == .true.) exit

      ! --- parse wall row: 6 reals + 1 int
      ios = 0
      read(line,*,iostat=ios) x1_left,x2_left,x3_left, x1_right,x2_right,x3_right, label_potential
      if (ios /= 0) error stop 'build_boundary: malformed wall row in geometry.inp'

      ! classify wall_type from signs BEFORE abs() (legacy semantics)
      if (label_potential < 0) then
        ! non-metal surface; periodic in Z if x3_right<0 (sign on x3_right)
        if (x3_right < 0.0_real64) then
          ! dielectric in YZ planes -> 3 or 4
          cnt_yz = cnt_yz + 1
          bnd%wall_type(abs(label_potential)) = merge(3,4, cnt_yz==1)
        else if (x2_right < 0.0_real64) then
          ! dielectric in XY planes -> code 2
          bnd%wall_type(abs(label_potential)) = 2
        else
          ! default non-metal tag (leave 0 => will still paint as -1)
          if (bnd%wall_type(abs(label_potential)) == 0) bnd%wall_type(abs(label_potential)) = 2
        end if
        if (x3_right < 0.0_real64) then
          bnd%flag_cylinder = .true.
          bnd%wall_type(abs(label_potential)) = -abs(bnd%wall_type(abs(label_potential)))  ! negative => “periodic with that surface”
        end if
      else
        ! metal
        if (bnd%wall_type(label_potential) == 0) bnd%wall_type(label_potential) = 1
      end if

      if (abs(label_potential) > bnd%nbPotentialLabel) error stop 'build_boundary: wall label > nbPotentialLabel'
      if (label_potential > 0) max_label = max(max_label, label_potential)

      ! clamp to box and convert cm -> m (legacy input)
      x1_left=max(0d0,min(x1_left,xmax*1d2))*1d-2; x1_right=max(0d0,min(x1_right,xmax*1d2))*1d-2
      x2_left=max(0d0,min(x2_left,ymax*1d2))*1d-2; x2_right=max(0d0,min(x2_right,ymax*1d2))*1d-2
      x3_left=max(0d0,min(x3_left,zmax*1d2))*1d-2; x3_right=max(0d0,min(x3_right,zmax*1d2))*1d-2

      ! indices on ghosted grid
      index_x1_left = int(x1_left/dx)+1; if (index_x1_left <= 1) index_x1_left = 0
      index_x1_right = int(x1_right/dx)+1; if (index_x1_right >= nb_cell_x1) index_x1_right = nb_cell_x1+2
      index_x2_left = int(x2_left/dy)+1; if (index_x2_left <= 1) index_x2_left = 0
      index_x2_right = int(x2_right/dy)+1; if (index_x2_right >= nb_cell_x2) index_x2_right = nb_cell_x2+2
      index_x3_left = int(x3_left/dz)+1; if (index_x3_left <= 1) index_x3_left = 0
      index_x3_right = int(x3_right/dz)+1; if (index_x3_right >= nb_cell_x3) index_x3_right = nb_cell_x3+2

      call paint_block(bnd, label_potential, index_x1_left,index_x1_right, index_x2_left,index_x2_right, index_x3_left,index_x3_right)
    end do
    close(iu)

    if (max_label > bnd%nbPotentialLabel) error stop 'build_boundary: label > nbPotentialLabel encountered'

    if (bnd%flag_cylinder) then
      ! clear ghost caps if periodic in Z
      bnd%bcnd(:,:,0:1)        = 0
      bnd%bcnd(:,:,nb_cell_x3+1:nb_cell_x3+2)  = 0
      bnd%phi (:,:,0:1)        = 0.0_real64
      bnd%phi (:,:,nb_cell_x3+1:nb_cell_x3+2)  = 0.0_real64
    end if

    ! ---- persist to HDF5 (state + metadata) ----
    origin  = [0.0_real64, 0.0_real64, 0.0_real64]
    spacing = [dx, dy, dz]
    nghost  = 2

    call save_array_h5("./Outputs/state.h5","bcnd", bnd%bcnd, origin, spacing, nghost, replace=.true.)
    call save_array_h5("./Outputs/state.h5","phi" , bnd%phi , origin, spacing, nghost, replace=.true.)

    ! metadata arrays (rank-1) – no origin/spacing needed but harmless if provided
    call save_array_h5("./Outputs/state.h5","wall_type", bnd%wall_type, replace=.true.)
    call save_array_h5("./Outputs/state.h5","potential_from_boundary"    , bnd%potential_from_boundary    , replace=.true.)

  end subroutine build_boundary

  !===================== helpers =====================

  subroutine paint_block(bnd, label_potential, index_x1_left,index_x1_right, index_x2_left,index_x2_right, index_x3_left,index_x3_right)
    type(Boundary), intent(inout) :: bnd
    integer, intent(in) :: label_potential, index_x1_left,index_x1_right, index_x2_left,index_x2_right, index_x3_left,index_x3_right
    integer :: ix,iy,iz, lab
    lab = abs(label_potential)

    do iz = max(0,index_x3_left), min(size(bnd%bcnd,3)-1,index_x3_right)
      do iy = max(0,index_x2_left), min(size(bnd%bcnd,2)-1,index_x2_right)
        do ix = max(0,index_x1_left), min(size(bnd%bcnd,1)-1,index_x1_right)
          if (label_potential >= 0) then                ! metal -> keep label number
            bnd%bcnd(ix,iy,iz) = lab
            if (lab > 0) bnd%phi(ix,iy,iz) = bnd%potential_from_boundary(lab)
          else                               ! any non-metal -> collapse to -1 (Option A)
            bnd%bcnd(ix,iy,iz) = -1
            bnd%phi (ix,iy,iz) = 0.0_real64
          end if
        end do
      end do
    end do
  end subroutine paint_block

end module mod_boundary
