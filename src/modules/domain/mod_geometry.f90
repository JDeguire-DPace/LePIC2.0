module mod_geometry
  use iso_fortran_env, only: real64, int16, int32
  use mpi
  use hdf5
  use mod_inputfiles, only: find_geometry_file
  use mod_functionsText
  implicit none
  private
  public :: Domain

  type :: Domain
    integer(int16) :: ncell_x1 = 1, ncell_x2 = 1, ncell_x3 = 1
    integer(int16) :: ncell_ghost = 2
    real(real64)   :: x1_min=0.0, x1_max=1.0, x2_min=0.0, x2_max=1.0, x3_min=0.0, x3_max=1.0
    real(real64)   :: Length_x1=1.0, Length_x2=1.0, Length_x3=1.0
    real(real64)   :: dx1=1.0, dx2=1.0, dx3=1.0
  contains
    procedure :: read_geometry    ! reads only line 1–2 (nx,ny,nz and extents)
    procedure :: get_nbPotentialLabel        ! reads last integer line = nbPotentialLabel
  end type Domain

  contains

  subroutine read_geometry(self)
    class(Domain), intent(inout) :: self
    character(len=256) :: gfile
    integer(int32) :: nx,ny,nz
    real(real64) :: xmax,ymax,zmax
    integer :: iu, ios

    gfile = trim(find_geometry_file())
    open(newunit=iu, file=gfile, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'read_geometry: cannot open '//gfile

    call next_data_line(iu, ios); if (ios /= 0) error stop 'read_geometry: missing nx ny nz'
    read(iu,*,iostat=ios) nx, ny, nz
    if (ios /= 0) error stop 'read_geometry: bad nx ny nz'

    call next_data_line(iu, ios); if (ios /= 0) error stop 'read_geometry: missing xmax ymax zmax'
    read(iu,*,iostat=ios) xmax, ymax, zmax
    if (ios /= 0) error stop 'read_geometry: bad xmax ymax zmax'
    close(iu)

    xmax = abs(xmax); ymax = abs(ymax); zmax = abs(zmax)
    ! legacy cm -> SI m
    xmax = xmax*1.0e-2; ymax = ymax*1.0e-2; zmax = zmax*1.0e-2

    self%ncell_x1 = int(nx,kind(self%ncell_x1))
    self%ncell_x2 = int(ny,kind(self%ncell_x2))
    self%ncell_x3 = int(nz,kind(self%ncell_x3))

    self%x1_min=0.0; self%x1_max=xmax
    self%x2_min=0.0; self%x2_max=ymax
    self%x3_min=0.0; self%x3_max=zmax

    self%Length_x1=xmax; self%Length_x2=ymax; self%Length_x3=zmax
    self%dx1 = self%Length_x1/real(self%ncell_x1,real64)
    self%dx2 = self%Length_x2/real(self%ncell_x2,real64)
    self%dx3 = self%Length_x3/real(self%ncell_x3,real64)
  end subroutine read_geometry

  integer function get_nbPotentialLabel(self) result(nbPotentialLabel)
    class(Domain), intent(in) :: self
    character(len=256) :: gfile, line
    integer :: iu, ios, tmp, max_label, ind
    real(real64) :: xl,yl,zl, xr,yr,zr

    gfile = trim(find_geometry_file())
    open(newunit=iu, file=gfile, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'get_nbPotentialLabel: cannot open '//gfile

    ! --- skip the first two data lines (nx,ny,nz) and (xmax,ymax,zmax) ---
    call next_data_line(iu, ios); if (ios /= 0) error stop 'get_nbPotentialLabel: missing header line 1'
    read(iu,*)
    call next_data_line(iu, ios); if (ios /= 0) error stop 'get_nbPotentialLabel: missing header line 2'
    read(iu,*)

    max_label = 0
    do
      call next_data_line(iu, ios); if (ios /= 0) exit
      read(iu,'(A)',iostat=ios) line
      if (ios /= 0) exit
      call strip_inline_comment(line)
      if (len_trim(line) == 0) cycle

      ! Try a full wall row: 6 reals + 1 integer label
      ios = 0
      read(line,*,iostat=ios) xl,yl,zl, xr,yr,zr, ind
      if (ios == 0) then
        max_label = max(max_label, abs(ind))
        cycle
      end if

      ! If it's a single integer (e.g. "5") at the end, IGNORE it and stop scanning.
      ios = 0
      read(line,*,iostat=ios) tmp
      if (ios == 0 .and. index(adjustl(line),' ') == 0) exit
      ! otherwise, keep scanning (malformed line or comment-only)
    end do
    close(iu)

    if (max_label <= 0) error stop 'get_nbPotentialLabel: could not infer nbPotentialLabel from wall rows'
    nbPotentialLabel = max_label
  end function get_nbPotentialLabel

end module mod_geometry
