module mod_MagneticField
  use iso_fortran_env, only: real64, int32, output_unit
  implicit none
  private

  !-----------------------------
  ! Public type & procedures
  !-----------------------------
  integer, parameter, public :: MFIELD_NONE     = 0
  integer, parameter, public :: MFIELD_GAUSSIAN = 1
  integer, parameter, public :: MFIELD_MAP      = 2

  type, public :: MagneticField
    ! --- grid ---
    integer(int32) :: nx1 = 1, nx2 = 1, nx3 = 1
    real(real64)   :: x1_min = 0.0d0, x1_max = 1.0d0
    real(real64)   :: x2_min = 0.0d0, x2_max = 1.0d0
    real(real64)   :: x3_min = 0.0d0, x3_max = 1.0d0

    ! --- model selection ---
    integer(int32) :: field_type  = MFIELD_NONE     ! gaussian|map
    integer(int32) :: direction   = 3               ! 1→x1, 2→x2, 3→x3

    ! --- gaussian parameters ---
    real(real64)   :: b0   = 0.0d0
    real(real64)   :: L    = 0.05d0
    real(real64)   :: x1_0 = 0.0d0, x2_0 = 0.0d0, x3_0 = 0.0d0

    ! --- external map parameters ---
    character(len=:), allocatable :: map_file
    real(real64)   :: map_scaling   = 1.0d0
    integer(int32) :: map_component = 3   ! which component the map adds to (default z)

    ! --- data ---
    real(real64), allocatable :: B(:,:,:,:)   ! (3, nx1, nx2, nx3)
  contains
    procedure, public :: init_defaults
    procedure, public :: read_from_file
    procedure, public :: allocate_arrays
    procedure, public :: build_field
    procedure, public :: grid_shape
    procedure, public :: spacing_and_origin
    procedure, public :: component => component_ptr
  end type MagneticField

  public :: decode_direction, decode_field_type

  contains
  !==============================
  subroutine init_defaults(self)
  !==============================
    class(MagneticField), intent(inout) :: self
    self%nx1 = 1; self%nx2 = 1; self%nx3 = 1
    self%x1_min = 0.0d0; self%x1_max = 1.0d0
    self%x2_min = 0.0d0; self%x2_max = 1.0d0
    self%x3_min = 0.0d0; self%x3_max = 1.0d0
    self%field_type = MFIELD_NONE
    self%direction  = 3
    self%b0   = 0.0d0
    self%L    = 0.05d0
    self%x1_0 = 0.0d0; self%x2_0 = 0.0d0; self%x3_0 = 0.0d0
    if (allocated(self%map_file)) deallocate(self%map_file)
    self%map_scaling   = 1.0d0
    self%map_component = 3
    if (allocated(self%B)) deallocate(self%B)
  end subroutine init_defaults

  !=========================================================
  subroutine read_from_file(self, filename)
  ! filename is OPTIONAL; falls back to find_magneticfield_file()
  !=========================================================
    use mod_InputFiles, only: find_magneticfield_file
    class(MagneticField), intent(inout) :: self
    character(len=*),     intent(in), optional :: filename

    character(len=:), allocatable :: path
    integer :: iu, ios
    character(len=:), allocatable :: line, key, val

    call self%init_defaults()

    ! Decide path
    if (present(filename)) then
      path = trim(filename)
    else
      path = trim(find_magneticfield_file())
    end if

    ! Open & parse
    open(newunit=iu, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(output_unit,'(a)') 'mod_MagneticField: cannot open '//path
      error stop 'mod_MagneticField: failed to open input'
    end if

    do
      call read_line_alloc(iu, line, ios); if (ios /= 0) exit
      call strip_inline_comment(line)
      if (is_blank(line)) cycle
      call parse_kv_eq(line, key, val)      ! key = value
      call to_lower_inplace(key)
      call trim_inplace(key); call trim_inplace(val)
      if (len_trim(key) == 0) cycle

      select case (key)
      case('nx1');           read(val,*) self%nx1
      case('nx2');           read(val,*) self%nx2
      case('nx3');           read(val,*) self%nx3
      case('x1_min');        read(val,*) self%x1_min
      case('x1_max');        read(val,*) self%x1_max
      case('x2_min');        read(val,*) self%x2_min
      case('x2_max');        read(val,*) self%x2_max
      case('x3_min');        read(val,*) self%x3_min
      case('x3_max');        read(val,*) self%x3_max
      case('direction');     self%direction  = decode_direction(val)
      case('field_type');    self%field_type = decode_field_type(val)
      case('b0');            read(val,*) self%b0
      case('l');             read(val,*) self%L
      case('x1_0');          read(val,*) self%x1_0
      case('x2_0');          read(val,*) self%x2_0
      case('x3_0');          read(val,*) self%x3_0
      case('map_file')
        if (allocated(self%map_file)) deallocate(self%map_file)
        allocate(character(len=len_trim(val)) :: self%map_file)
        self%map_file = trim(val)
      case('map_scaling');   read(val,*) self%map_scaling
      case('map_component'); self%map_component = decode_direction(val)
      case default
        ! ignore unknown keys
      end select
    end do
    close(iu)
  end subroutine read_from_file

  !==================================
  subroutine allocate_arrays(self)
  !==================================
    class(MagneticField), intent(inout) :: self
    if (allocated(self%B)) deallocate(self%B)
    allocate(self%B(3, self%nx1, self%nx2, self%nx3))
    self%B = 0.0d0
  end subroutine allocate_arrays

  !=============================
  subroutine build_field(self)
  !=============================
    class(MagneticField), intent(inout) :: self
    integer :: i, j, k, comp
    real(real64) :: dx1, dx2, dx3, x, y, z

    if (.not.allocated(self%B)) call self%allocate_arrays()
    self%B = 0.0d0

    dx1 = merge( (self%x1_max-self%x1_min)/real(max(1,self%nx1-1),real64), 0.0d0, self%nx1>1 )
    dx2 = merge( (self%x2_max-self%x2_min)/real(max(1,self%nx2-1),real64), 0.0d0, self%nx2>1 )
    dx3 = merge( (self%x3_max-self%x3_min)/real(max(1,self%nx3-1),real64), 0.0d0, self%nx3>1 )

    select case (self%field_type)
    case (MFIELD_GAUSSIAN)
      comp = self%direction
      do k = 1, self%nx3
        z = self%x3_min + dx3*real(k-1,real64)
        do j = 1, self%nx2
          y = self%x2_min + dx2*real(j-1,real64)
          do i = 1, self%nx1
            x = self%x1_min + dx1*real(i-1,real64)
            self%B(comp,i,j,k) = self%b0 * exp( -((x-self%x1_0)**2 + (y-self%x2_0)**2 + (z-self%x3_0)**2) / (2.0d0*self%L**2) )
          end do
        end do
      end do

    case (MFIELD_MAP)
      call apply_map_to_component(self)

    case default
      ! nothing
    end select
  end subroutine build_field

  !-------------------------------------
  subroutine apply_map_to_component(self)
  !-------------------------------------
    class(MagneticField), intent(inout) :: self
    integer :: iu, ios
    integer :: ix, iy, iz
    real(real64) :: bval
    integer :: comp
    character(len=:), allocatable :: path

    if (.not.allocated(self%map_file)) then
      write(output_unit,'(a)') 'mod_MagneticField: map field selected but no map_file provided'
      error stop 'mod_MagneticField: missing map_file'
    end if
    path = trim(self%map_file)
    comp = self%map_component

    open(newunit=iu, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(output_unit,'(a)') 'mod_MagneticField: cannot open map file '//path
      error stop 'mod_MagneticField: cannot open map file'
    end if

    do
      read(iu,*,iostat=ios) ix, iy, iz, bval
      if (ios /= 0) exit
      if (ix>=1 .and. ix<=self%nx1 .and. iy>=1 .and. iy<=self%nx2 .and. iz>=1 .and. iz<=self%nx3) then
        self%B(comp,ix,iy,iz) = self%B(comp,ix,iy,iz) + self%map_scaling*bval
      end if
    end do
    close(iu)
  end subroutine apply_map_to_component

  !=============================
  subroutine grid_shape(self, n1,n2,n3)
  !=============================
    class(MagneticField), intent(in) :: self
    integer(int32), intent(out) :: n1,n2,n3
    n1 = self%nx1; n2 = self%nx2; n3 = self%nx3
  end subroutine grid_shape

  !===============================================
  subroutine spacing_and_origin(self, spacing, origin)
  !===============================================
    class(MagneticField), intent(in) :: self
    real(real64), intent(out) :: spacing(3), origin(3)
    spacing(1) = merge( (self%x1_max-self%x1_min)/real(max(1,self%nx1-1),real64), 1.0d0, self%nx1>1 )
    spacing(2) = merge( (self%x2_max-self%x2_min)/real(max(1,self%nx2-1),real64), 1.0d0, self%nx2>1 )
    spacing(3) = merge( (self%x3_max-self%x3_min)/real(max(1,self%nx3-1),real64), 1.0d0, self%nx3>1 )
    origin(1)  = self%x1_min
    origin(2)  = self%x2_min
    origin(3)  = self%x3_min
  end subroutine spacing_and_origin

  !========================================
  function component_ptr(self, comp) result(view)
  ! Return a pointer to B(comp,:,:,:)
  !========================================
    class(MagneticField), target, intent(inout) :: self
    integer, intent(in) :: comp
    real(real64), pointer :: view(:,:,:)
    if (.not.allocated(self%B)) nullify(view)
    if (allocated(self%B)) then
      view => self%B(comp, :, :, :)
    end if
  end function component_ptr

  !==========================
  pure integer function decode_direction(s)
  !==========================
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: t
    t = trim(adjustl(s)); call to_lower_inplace(t)
    select case (t)
    case('1','x','x1','b1'); decode_direction = 1
    case('2','y','x2','b2'); decode_direction = 2
    case('3','z','x3','b3'); decode_direction = 3
    case default;             decode_direction = 3
    end select
  end function decode_direction

  !==========================
  pure integer function decode_field_type(s)
  !==========================
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: t
    t = trim(adjustl(s)); call to_lower_inplace(t)
    select case (t)
    case('gaussian'); decode_field_type = MFIELD_GAUSSIAN
    case('map');      decode_field_type = MFIELD_MAP
    case default;     decode_field_type = MFIELD_NONE
    end select
  end function decode_field_type

  !================ Helpers (local, minimal) ================

  subroutine read_line_alloc(iu, line, ios)
    integer, intent(in) :: iu
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out) :: ios
    character(len=1024) :: tmp
    read(iu,'(A)',iostat=ios) tmp
    if (ios == 0) then
      line = trim(tmp)
    else
      if (allocated(line)) deallocate(line)
    end if
  end subroutine read_line_alloc

  pure subroutine strip_inline_comment(s)
    character(len=*), intent(inout) :: s
    integer :: p
    p = index(s, '#'); if (p==0) p = index(s, '!')
    if (p > 0) s = s(:p-1)
  end subroutine strip_inline_comment

  pure logical function is_blank(s)
    character(len=*), intent(in) :: s
    is_blank = (len_trim(adjustl(s)) == 0)
  end function is_blank

  pure subroutine parse_kv_eq(line, key, val)
    character(len=*), intent(in)  :: line
    character(len=:), allocatable, intent(out) :: key, val
    integer :: p
    character(len=:), allocatable :: L
    L = trim(line)
    p = index(L, '=')
    if (p <= 0) then
      key = ''; val = ''
    else
      key = trim(L(:p-1))
      val = trim(L(p+1:))
    end if
  end subroutine parse_kv_eq

  pure subroutine to_lower_inplace(s)
    character(len=*), intent(inout) :: s
    integer :: i, ia
    do i=1,len(s)
      ia = iachar(s(i:i))
      if (ia>=iachar('A') .and. ia<=iachar('Z')) s(i:i) = achar(ia+32)
    end do
  end subroutine to_lower_inplace

  pure subroutine trim_inplace(s)
    character(len=*), intent(inout) :: s
    s = trim(adjustl(s))
  end subroutine trim_inplace

end module mod_MagneticField
