module mod_io_hdf5
  use iso_fortran_env, only: int32, real64
  use hdf5
  implicit none
  private

  ! Public generics
  public :: save_array_h5, load_array_h5
  public :: save_bcnd_h5,  load_bcnd_h5

  ! ---------- Generic bindings (must be at module scope, before CONTAINS) ----------
  interface save_array_h5
     module procedure save_i32_1d, save_i32_2d, save_i32_3d
     module procedure save_r64_1d, save_r64_2d, save_r64_3d
  end interface

  interface load_array_h5
     module procedure load_i32_1d, load_i32_2d, load_i32_3d
     module procedure load_r64_1d, load_r64_2d, load_r64_3d
  end interface
  ! -------------------------------------------------------------------------------

  contains

  ! ==========================
  ! === Core save helpers  ===
  ! ==========================
  subroutine save_core_int(filename, dset, A, dims_h5, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    integer(int32),             intent(in) :: A(..)         ! rank-agnostic actuals will pass via wrappers
    integer(HSIZE_T),           intent(in) :: dims_h5(:)    ! size = rank
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace

    integer(HID_T) :: file_id, dset_id, dspace_id, attr_id
    integer        :: ierr, rank_h5, ngh, tcode
    logical :: file_opened, exists, want_replace

    rank_h5 = size(dims_h5)
    ngh = 0; if (present(nghost)) ngh = nghost
    want_replace = .false.; if (present(replace)) want_replace = replace
    tcode = 0  ! 0=int32, 1=real64

    call h5open_f(ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr)
    file_opened = (ierr == 0)
    if (.not. file_opened) then
      call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr)
      if (ierr /= 0) error stop 'save_array_h5: h5fcreate_f failed'
    end if

    call h5lexists_f(file_id, trim(dset), exists, ierr)
    if (exists) then
      if (want_replace) then
        call h5ldelete_f(file_id, trim(dset), ierr)
        if (ierr /= 0) error stop 'save_array_h5: failed to delete existing dataset'
      else
        error stop 'save_array_h5: dataset exists; pass replace=.true. to overwrite'
      end if
    end if

    call h5screate_simple_f(rank_h5, dims_h5, dspace_id, ierr)
    if (ierr /= 0) error stop 'save_array_h5: h5screate_simple_f failed'

    call h5dcreate_f(file_id, trim(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
    if (ierr /= 0) error stop 'save_array_h5: h5dcreate_f failed (int)'

    select rank (A)
    rank (1); call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A, dims_h5, ierr)
    rank (2); call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A, dims_h5, ierr)
    rank (3); call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A, dims_h5, ierr)
    rank default
      error stop 'save_array_h5(int): unsupported rank'
    end select
    if (ierr /= 0) error stop 'save_array_h5: h5dwrite_f failed'

    call write_attrs(dset_id, rank_h5, dims_h5, tcode, origin, spacing, ngh)

    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
  end subroutine save_core_int


  subroutine save_core_real(filename, dset, A, dims_h5, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(..)
    integer(HSIZE_T),           intent(in) :: dims_h5(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace

    integer(HID_T) :: file_id, dset_id, dspace_id, attr_id
    integer        :: ierr, rank_h5, ngh, tcode
    logical :: file_opened, exists, want_replace

    rank_h5 = size(dims_h5)
    ngh = 0; if (present(nghost)) ngh = nghost
    want_replace = .false.; if (present(replace)) want_replace = replace
    tcode = 1  ! 0=int32, 1=real64

    call h5open_f(ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr)
    file_opened = (ierr == 0)
    if (.not. file_opened) then
      call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr)
      if (ierr /= 0) error stop 'save_array_h5: h5fcreate_f failed'
    end if

    call h5lexists_f(file_id, trim(dset), exists, ierr)
    if (exists) then
      if (want_replace) then
        call h5ldelete_f(file_id, trim(dset), ierr)
        if (ierr /= 0) error stop 'save_array_h5: failed to delete existing dataset'
      else
        error stop 'save_array_h5: dataset exists; pass replace=.true. to overwrite'
      end if
    end if

    call h5screate_simple_f(rank_h5, dims_h5, dspace_id, ierr)
    if (ierr /= 0) error stop 'save_array_h5: h5screate_simple_f failed'

    call h5dcreate_f(file_id, trim(dset), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
    if (ierr /= 0) error stop 'save_array_h5: h5dcreate_f failed (real)'

    select rank (A)
    rank (1); call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, dims_h5, ierr)
    rank (2); call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, dims_h5, ierr)
    rank (3); call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, dims_h5, ierr)
    rank default
      error stop 'save_array_h5(real): unsupported rank'
    end select
    if (ierr /= 0) error stop 'save_array_h5: h5dwrite_f failed'

    call write_attrs(dset_id, rank_h5, dims_h5, tcode, origin, spacing, ngh)

    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
  end subroutine save_core_real


  subroutine write_attrs(dset_id, rank_h5, dims_h5, dtype_code, origin, spacing, nghost)
    integer(HID_T),             intent(in) :: dset_id
    integer,                    intent(in) :: rank_h5, dtype_code
    integer(HSIZE_T),           intent(in) :: dims_h5(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost

    integer(HID_T) :: as1, as3, attr_id
    integer        :: ierr
    integer(HSIZE_T) :: one(1), three(1)
    integer(int32) :: v3(3)

    one  = [1_HSIZE_T]
    three= [3_HSIZE_T]

    ! len-3 dataspace
    call h5screate_simple_f(1, three, as3, ierr)
    if (present(origin)) then
      call h5acreate_f(dset_id, "origin",  H5T_NATIVE_DOUBLE, as3, attr_id, ierr)
      call h5awrite_f (attr_id, H5T_NATIVE_DOUBLE, origin, three, ierr)
      call h5aclose_f(attr_id, ierr)
    end if
    if (present(spacing)) then
      call h5acreate_f(dset_id, "spacing", H5T_NATIVE_DOUBLE, as3, attr_id, ierr)
      call h5awrite_f (attr_id, H5T_NATIVE_DOUBLE, spacing, three, ierr)
      call h5aclose_f(attr_id, ierr)
    end if

    v3 = 1
    if (rank_h5 >= 1) v3(1) = int(dims_h5(1))
    if (rank_h5 >= 2) v3(2) = int(dims_h5(2))
    if (rank_h5 >= 3) v3(3) = int(dims_h5(3))
    call h5acreate_f(dset_id, "nxyz", H5T_NATIVE_INTEGER, as3, attr_id, ierr)
    call h5awrite_f (attr_id, H5T_NATIVE_INTEGER, v3, three, ierr)
    call h5aclose_f(attr_id, ierr)

    v3 = [1,0,0]
    call h5acreate_f(dset_id, "format_version_semver", H5T_NATIVE_INTEGER, as3, attr_id, ierr)
    call h5awrite_f (attr_id, H5T_NATIVE_INTEGER, v3, three, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5sclose_f(as3, ierr)

    ! len-1 dataspace
    call h5screate_simple_f(1, one, as1, ierr)
    call h5acreate_f(dset_id, "dtype", H5T_NATIVE_INTEGER, as1, attr_id, ierr)
    call h5awrite_f (attr_id, H5T_NATIVE_INTEGER, dtype_code, one, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5acreate_f(dset_id, "rank",  H5T_NATIVE_INTEGER, as1, attr_id, ierr)
    call h5awrite_f (attr_id, H5T_NATIVE_INTEGER, rank_h5,   one, ierr)
    call h5aclose_f(attr_id, ierr)

    if (present(nghost)) then
      call h5acreate_f(dset_id, "nghost", H5T_NATIVE_INTEGER, as1, attr_id, ierr)
      call h5awrite_f (attr_id, H5T_NATIVE_INTEGER, nghost,   one, ierr)
      call h5aclose_f(attr_id, ierr)
    end if

    call h5sclose_f(as1, ierr)
  end subroutine write_attrs


  ! ==========================
  ! === Typed save wrappers ==
  ! ==========================
  subroutine save_i32_1d(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    integer(int32),             intent(in) :: A(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    call save_core_int(filename, dset, A, [ int(size(A,1),HSIZE_T) ], origin, spacing, nghost, replace)
  end subroutine

  subroutine save_i32_2d(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    integer(int32),             intent(in) :: A(:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    call save_core_int(filename, dset, A, [ int(size(A,1),HSIZE_T), int(size(A,2),HSIZE_T) ], origin, spacing, nghost, replace)
  end subroutine

  subroutine save_i32_3d(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    integer(int32),             intent(in) :: A(:,:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    call save_core_int(filename, dset, A, [ int(size(A,1),HSIZE_T), int(size(A,2),HSIZE_T), int(size(A,3),HSIZE_T) ], origin, spacing, nghost, replace)
  end subroutine

  subroutine save_r64_1d(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    call save_core_real(filename, dset, A, [ int(size(A,1),HSIZE_T) ], origin, spacing, nghost, replace)
  end subroutine

  subroutine save_r64_2d(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    call save_core_real(filename, dset, A, [ int(size(A,1),HSIZE_T), int(size(A,2),HSIZE_T) ], origin, spacing, nghost, replace)
  end subroutine

  subroutine save_r64_3d(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(:,:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    call save_core_real(filename, dset, A, [ int(size(A,1),HSIZE_T), int(size(A,2),HSIZE_T), int(size(A,3),HSIZE_T) ], origin, spacing, nghost, replace)
  end subroutine


  ! ==========================
  ! === Core load helpers  ===
  ! ==========================
  subroutine load_core_int(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    integer(int32), allocatable, intent(out) :: A(..)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost

    integer(HID_T) :: file_id, dset_id, dspace_id, attr_id
    integer        :: ierr, rank
    integer(HSIZE_T), allocatable :: dims(:), maxdims(:)

    call h5open_f(ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr); if (ierr/=0) error stop 'load_array_h5: h5fopen_f failed'
    call h5dopen_f(file_id, trim(dset), dset_id, ierr);              if (ierr/=0) error stop 'load_array_h5: h5dopen_f failed'
    call h5dget_space_f(dset_id, dspace_id, ierr)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
    allocate(dims(rank), maxdims(rank))
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)

    select rank (A)
    rank (1)
      allocate(A( int(dims(1)) ))
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, A, dims, ierr)
    rank (2)
      allocate(A( int(dims(1)), int(dims(2)) ))
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, A, dims, ierr)
    rank (3)
      allocate(A( int(dims(1)), int(dims(2)), int(dims(3)) ))
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, A, dims, ierr)
    rank default
      error stop 'load_array_h5(int32): unexpected rank'
    end select

    call read_common_attrs(dset_id, origin, spacing, nghost)
    call h5sclose_f(dspace_id, ierr); call h5dclose_f(dset_id, ierr); call h5fclose_f(file_id, ierr); call h5close_f(ierr)
  end subroutine load_core_int


  subroutine load_core_real(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    real(real64),   allocatable, intent(out) :: A(..)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost

    integer(HID_T) :: file_id, dset_id, dspace_id, attr_id
    integer        :: ierr, rank
    integer(HSIZE_T), allocatable :: dims(:), maxdims(:)

    call h5open_f(ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr); if (ierr/=0) error stop 'load_array_h5: h5fopen_f failed'
    call h5dopen_f(file_id, trim(dset), dset_id, ierr);              if (ierr/=0) error stop 'load_array_h5: h5dopen_f failed'
    call h5dget_space_f(dset_id, dspace_id, ierr)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
    allocate(dims(rank), maxdims(rank))
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)

    select rank (A)
    rank (1)
      allocate(A( int(dims(1)) ))
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, ierr)
    rank (2)
      allocate(A( int(dims(1)), int(dims(2)) ))
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, ierr)
    rank (3)
      allocate(A( int(dims(1)), int(dims(2)), int(dims(3)) ))
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, ierr)
    rank default
      error stop 'load_array_h5(real64): unexpected rank'
    end select

    call read_common_attrs(dset_id, origin, spacing, nghost)
    call h5sclose_f(dspace_id, ierr); call h5dclose_f(dset_id, ierr); call h5fclose_f(file_id, ierr); call h5close_f(ierr)
  end subroutine load_core_real


  subroutine read_common_attrs(dset_id, origin, spacing, nghost)
    integer(HID_T),             intent(in)  :: dset_id
    real(real64), dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost
    integer(HID_T) :: attr_id
    integer        :: ierr
    logical :: has
    integer(HSIZE_T) :: one(1), three(1)
    one   = [1_HSIZE_T]
    three = [3_HSIZE_T]

    if (present(origin)) then
      call h5aexists_f(dset_id, "origin", has, ierr)
      if (has) then
        call h5aopen_f(dset_id, "origin", attr_id, ierr)
        call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, origin, three, ierr)
        call h5aclose_f(attr_id, ierr)
      else
        origin = 0.0_real64
      end if
    end if
    if (present(spacing)) then
      call h5aexists_f(dset_id, "spacing", has, ierr)
      if (has) then
        call h5aopen_f(dset_id, "spacing", attr_id, ierr)
        call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, spacing, three, ierr)
        call h5aclose_f(attr_id, ierr)
      else
        spacing = 0.0_real64
      end if
    end if
    if (present(nghost)) then
      call h5aexists_f(dset_id, "nghost", has, ierr)
      if (has) then
        call h5aopen_f(dset_id, "nghost", attr_id, ierr)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, nghost, one, ierr)
        call h5aclose_f(attr_id, ierr)
      else
        nghost = 0
      end if
    end if
  end subroutine read_common_attrs


  ! ==========================
  ! === Typed load wrappers ==
  ! ==========================
  subroutine load_i32_1d(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    integer(int32), allocatable, intent(out) :: A(:)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost
    call load_core_int(filename, dset, A, origin, spacing, nghost)
  end subroutine

  subroutine load_i32_2d(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    integer(int32), allocatable, intent(out) :: A(:,:)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost
    call load_core_int(filename, dset, A, origin, spacing, nghost)
  end subroutine

  subroutine load_i32_3d(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    integer(int32), allocatable, intent(out) :: A(:,:,:)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost
    call load_core_int(filename, dset, A, origin, spacing, nghost)
  end subroutine

  subroutine load_r64_1d(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    real(real64),   allocatable, intent(out) :: A(:)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost
    call load_core_real(filename, dset, A, origin, spacing, nghost)
  end subroutine

  subroutine load_r64_2d(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    real(real64),   allocatable, intent(out) :: A(:,:)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost
    call load_core_real(filename, dset, A, origin, spacing, nghost)
  end subroutine

  subroutine load_r64_3d(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    real(real64),   allocatable, intent(out) :: A(:,:,:)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost
    call load_core_real(filename, dset, A, origin, spacing, nghost)
  end subroutine


  ! ==========================
  ! === Back-compat wrappers =
  ! ==========================
  subroutine save_bcnd_h5(filename, bcnd, origin, spacing, nghost)
    character(*),               intent(in) :: filename
    integer(int32),             intent(in) :: bcnd(:,:,:)
    real(real64), dimension(3), intent(in) :: origin, spacing
    integer(int32),             intent(in) :: nghost
    call save_array_h5(filename, "bcnd", bcnd, origin, spacing, nghost, replace=.true.)
  end subroutine save_bcnd_h5

  subroutine load_bcnd_h5(filename, bcnd, origin, spacing, nghost)
    character(*),                intent(in)  :: filename
    integer(int32), allocatable, intent(out) :: bcnd(:,:,:)
    real(real64),   dimension(3), intent(out) :: origin, spacing
    integer(int32),               intent(out) :: nghost
    call load_array_h5(filename, "bcnd", bcnd, origin, spacing, nghost)
  end subroutine load_bcnd_h5

end module mod_io_hdf5
