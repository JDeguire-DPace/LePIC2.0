module mod_io_hdf5
  use iso_fortran_env, only: int32, real64
  use mpi
  use hdf5
  implicit none
  private

  ! Public generics
  public :: save_array_h5, load_array_h5
  public :: save_bcnd_h5,  load_bcnd_h5
  public :: save_array_h5_wo_ghost

  ! ---------- Generic bindings (must be at module scope, before CONTAINS) ----------
  interface save_array_h5
     module procedure save_i32_1d, save_i32_2d, save_i32_3d
     module procedure save_r64_1d, save_r64_2d, save_r64_3d
  end interface

  interface save_array_h5_wo_ghost
     module procedure save_i32_1d_wo_ghost, save_i32_2d_wo_ghost, save_i32_3d_wo_ghost
     module procedure save_r64_1d_wo_ghost, save_r64_2d_wo_ghost, save_r64_3d_wo_ghost
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
    integer(int32),             intent(in) :: A(..)
    integer(HSIZE_T),           intent(in) :: dims_h5(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace

    integer(HID_T) :: file_id, dset_id, dspace_id
    integer        :: ierr, rank_h5, ngh, tcode
    logical :: exists, want_replace, dset_exists
    integer        :: myrank, mpier

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpier)

    rank_h5 = size(dims_h5)
    ngh = 0; if (present(nghost)) ngh = nghost
    want_replace = .false.; if (present(replace)) want_replace = replace
    tcode = 0  ! 0=int32, 1=real64

    if (myrank == 0) then
      call h5open_f(ierr)

      ! Open or create file on rank 0 only
      inquire(file=trim(filename), exist=exists)
      if (exists) then
        call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr)
        if (ierr /= 0) error stop 'save_array_h5(int): h5fopen_f failed'
      else
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr)
        if (ierr /= 0) error stop 'save_array_h5(int): h5fcreate_f failed'
      end if

      ! delete if present and replace requested
      call h5lexists_f(file_id, trim(dset), dset_exists, ierr)
      if (dset_exists) then
        if (want_replace) then
            call h5ldelete_f(file_id, trim(dset), ierr)
            if (ierr /= 0) error stop 'save_array_h5(int): h5ldelete_f failed'
        else
            error stop 'save_array_h5(int): dataset exists; set replace=.true.'
        end if
      end if

      call h5screate_simple_f(rank_h5, dims_h5, dspace_id, ierr)
      if (ierr /= 0) error stop 'save_array_h5(int): h5screate_simple_f failed'

      call h5dcreate_f(file_id, trim(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
      if (ierr /= 0) error stop 'save_array_h5(int): h5dcreate_f failed'

      select rank (A)
      rank (1); call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A, dims_h5, ierr)
      rank (2); call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A, dims_h5, ierr)
      rank (3); call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A, dims_h5, ierr)
      rank default
        error stop 'save_array_h5(int): unsupported rank'
      end select
      if (ierr /= 0) error stop 'save_array_h5(int): h5dwrite_f failed'

      call write_attrs(dset_id, rank_h5, dims_h5, tcode, origin, spacing, ngh)

      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
    end if

    call MPI_Barrier(MPI_COMM_WORLD, mpier)
  end subroutine save_core_int



  subroutine save_core_real(filename, dset, A, dims_h5, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(..)
    integer(HSIZE_T),           intent(in) :: dims_h5(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace

    integer(HID_T) :: file_id, dset_id, dspace_id
    integer        :: ierr, rank_h5, ngh, tcode
    logical :: exists, want_replace, dset_exists
    integer        :: myrank, mpier

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpier)

    rank_h5 = size(dims_h5)
    ngh = 0; if (present(nghost)) ngh = nghost
    want_replace = .false.; if (present(replace)) want_replace = replace
    tcode = 1  ! 0=int32, 1=real64

    if (myrank == 0) then
      call h5open_f(ierr)

      inquire(file=trim(filename), exist=exists)
      if (exists) then
        call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, ierr)
        if (ierr /= 0) error stop 'save_array_h5(real): h5fopen_f failed'
      else
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr)
        if (ierr /= 0) error stop 'save_array_h5(real): h5fcreate_f failed'
      end if

      call h5lexists_f(file_id, trim(dset), dset_exists, ierr)
      if (dset_exists) then
        if (want_replace) then
            call h5ldelete_f(file_id, trim(dset), ierr)
            if (ierr /= 0) error stop 'save_array_h5(real): h5ldelete_f failed'
        else
            error stop 'save_array_h5(real): dataset exists; set replace=.true.'
        end if
      end if

      call h5screate_simple_f(rank_h5, dims_h5, dspace_id, ierr)
      if (ierr /= 0) error stop 'save_array_h5(real): h5screate_simple_f failed'

      call h5dcreate_f(file_id, trim(dset), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
      if (ierr /= 0) error stop 'save_array_h5(real): h5dcreate_f failed'

      select rank (A)
      rank (1); call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, dims_h5, ierr)
      rank (2); call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, dims_h5, ierr)
      rank (3); call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, dims_h5, ierr)
      rank default
        error stop 'save_array_h5(real): unsupported rank'
      end select
      if (ierr /= 0) error stop 'save_array_h5(real): h5dwrite_f failed'

      call write_attrs(dset_id, rank_h5, dims_h5, tcode, origin, spacing, ngh)

      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
    end if

    call MPI_Barrier(MPI_COMM_WORLD, mpier)
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

  ! =========================================================
  ! === Save WITHOUT ghost cells ===========================
  ! =========================================================

  subroutine save_i32_1d_wo_ghost(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    integer(int32),             intent(in) :: A(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    integer           :: ngh, iL, iU, n_tot
    integer(HSIZE_T)  :: dims_h5(1)
    integer(int32), allocatable :: tmp(:)

    ngh = 0; if (present(nghost)) ngh = nghost
    n_tot = ubound(A,1) - lbound(A,1) + 1
    if (2*ngh >= n_tot) error stop 'save_i32_1d_wo_ghost: nghost too large'

    iL = lbound(A,1) + ngh
    iU = ubound(A,1) - ngh
    allocate(tmp(iU-iL+1))
    tmp = A(iL:iU)
    dims_h5(1) = int(size(tmp,1), HSIZE_T)

    if (present(replace)) then
      call save_core_int(filename, dset, tmp, dims_h5, origin, spacing, 0_int32, replace)
    else
      call save_core_int(filename, dset, tmp, dims_h5, origin, spacing, 0_int32)
    end if
    deallocate(tmp)
  end subroutine save_i32_1d_wo_ghost


  subroutine save_i32_2d_wo_ghost(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    integer(int32),             intent(in) :: A(:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    integer           :: ngh
    integer           :: iL, iU, jL, jU, nx_tot, ny_tot
    integer(HSIZE_T)  :: dims_h5(2)
    integer(int32), allocatable :: tmp(:,:)

    ngh = 0; if (present(nghost)) ngh = nghost
    nx_tot = ubound(A,1) - lbound(A,1) + 1
    ny_tot = ubound(A,2) - lbound(A,2) + 1
    if (2*ngh >= nx_tot .or. 2*ngh >= ny_tot) error stop 'save_i32_2d_wo_ghost: nghost too large'

    iL = lbound(A,1) + ngh; iU = ubound(A,1) - ngh
    jL = lbound(A,2) + ngh; jU = ubound(A,2) - ngh
    allocate(tmp(iU-iL+1, jU-jL+1))
    tmp = A(iL:iU, jL:jU)
    dims_h5 = [int(size(tmp,1),HSIZE_T), int(size(tmp,2),HSIZE_T)]

    if (present(replace)) then
      call save_core_int(filename, dset, tmp, dims_h5, origin, spacing, 0_int32, replace)
    else
      call save_core_int(filename, dset, tmp, dims_h5, origin, spacing, 0_int32)
    end if
    deallocate(tmp)
  end subroutine save_i32_2d_wo_ghost


  subroutine save_i32_3d_wo_ghost(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    integer(int32),             intent(in) :: A(:,:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    integer           :: ngh, iL, iU, jL, jU, kL, kU
    integer           :: nx_tot, ny_tot, nz_tot
    integer(HSIZE_T)  :: dims_h5(3)
    integer(int32), allocatable :: tmp(:,:,:)

    ngh = 0; if (present(nghost)) ngh = nghost
    nx_tot = ubound(A,1)-lbound(A,1)+1; ny_tot = ubound(A,2)-lbound(A,2)+1; nz_tot = ubound(A,3)-lbound(A,3)+1
    if (2*ngh >= nx_tot .or. 2*ngh >= ny_tot .or. 2*ngh >= nz_tot) error stop 'save_i32_3d_wo_ghost: nghost too large'

    iL=lbound(A,1)+ngh; iU=ubound(A,1)-ngh
    jL=lbound(A,2)+ngh; jU=ubound(A,2)-ngh
    kL=lbound(A,3)+ngh; kU=ubound(A,3)-ngh

    allocate(tmp(iU-iL+1,jU-jL+1,kU-kL+1))
    tmp = A(iL:iU, jL:jU, kL:kU)
    dims_h5 = [int(size(tmp,1),HSIZE_T), int(size(tmp,2),HSIZE_T), int(size(tmp,3),HSIZE_T)]

    if (present(replace)) then
      call save_core_int(filename, dset, tmp, dims_h5, origin, spacing, 0_int32, replace)
    else
      call save_core_int(filename, dset, tmp, dims_h5, origin, spacing, 0_int32)
    end if
    deallocate(tmp)
  end subroutine save_i32_3d_wo_ghost


  subroutine save_r64_1d_wo_ghost(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    integer           :: ngh, iL, iU, n_tot
    integer(HSIZE_T)  :: dims_h5(1)
    real(real64), allocatable :: tmp(:)

    ngh = 0; if (present(nghost)) ngh = nghost
    n_tot = ubound(A,1) - lbound(A,1) + 1
    if (2*ngh >= n_tot) error stop 'save_r64_1d_wo_ghost: nghost too large'

    iL = lbound(A,1) + ngh; iU = ubound(A,1) - ngh
    allocate(tmp(iU-iL+1))
    tmp = A(iL:iU)
    dims_h5(1) = int(size(tmp,1), HSIZE_T)

    if (present(replace)) then
      call save_core_real(filename, dset, tmp, dims_h5, origin, spacing, 0_int32, replace)
    else
      call save_core_real(filename, dset, tmp, dims_h5, origin, spacing, 0_int32)
    end if
    deallocate(tmp)
  end subroutine save_r64_1d_wo_ghost


  subroutine save_r64_2d_wo_ghost(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    integer           :: ngh, iL, iU, jL, jU, nx_tot, ny_tot
    integer(HSIZE_T)  :: dims_h5(2)
    real(real64), allocatable :: tmp(:,:)

    ngh = 0; if (present(nghost)) ngh = nghost
    nx_tot = ubound(A,1)-lbound(A,1)+1; ny_tot = ubound(A,2)-lbound(A,2)+1
    if (2*ngh >= nx_tot .or. 2*ngh >= ny_tot) error stop 'save_r64_2d_wo_ghost: nghost too large'

    iL=lbound(A,1)+ngh; iU=ubound(A,1)-ngh
    jL=lbound(A,2)+ngh; jU=ubound(A,2)-ngh
    allocate(tmp(iU-iL+1,jU-jL+1))
    tmp = A(iL:iU, jL:jU)
    dims_h5 = [int(size(tmp,1),HSIZE_T), int(size(tmp,2),HSIZE_T)]

    if (present(replace)) then
      call save_core_real(filename, dset, tmp, dims_h5, origin, spacing, 0_int32, replace)
    else
      call save_core_real(filename, dset, tmp, dims_h5, origin, spacing, 0_int32)
    end if
    deallocate(tmp)
  end subroutine save_r64_2d_wo_ghost


  subroutine save_r64_3d_wo_ghost(filename, dset, A, origin, spacing, nghost, replace)
    character(*),               intent(in) :: filename, dset
    real(real64),               intent(in) :: A(:,:,:)
    real(real64), dimension(3), intent(in), optional :: origin, spacing
    integer(int32),             intent(in), optional :: nghost
    logical,                    intent(in), optional :: replace
    integer           :: ngh, iL, iU, jL, jU, kL, kU
    integer           :: nx_tot, ny_tot, nz_tot
    integer(HSIZE_T)  :: dims_h5(3)
    real(real64), allocatable :: tmp(:,:,:)

    ngh = 0; if (present(nghost)) ngh = nghost
    nx_tot = ubound(A,1)-lbound(A,1)+1; ny_tot = ubound(A,2)-lbound(A,2)+1; nz_tot = ubound(A,3)-lbound(A,3)+1
    if (2*ngh >= nx_tot .or. 2*ngh >= ny_tot .or. 2*ngh >= nz_tot) error stop 'save_r64_3d_wo_ghost: nghost too large'

    iL=lbound(A,1)+ngh; iU=ubound(A,1)-ngh
    jL=lbound(A,2)+ngh; jU=ubound(A,2)-ngh
    kL=lbound(A,3)+ngh; kU=ubound(A,3)-ngh

    allocate(tmp(iU-iL+1,jU-jL+1,kU-kL+1))
    tmp = A(iL:iU, jL:jU, kL:kU)
    dims_h5 = [int(size(tmp,1),HSIZE_T), int(size(tmp,2),HSIZE_T), int(size(tmp,3),HSIZE_T)]

    if (present(replace)) then
      call save_core_real(filename, dset, tmp, dims_h5, origin, spacing, 0_int32, replace)
    else
      call save_core_real(filename, dset, tmp, dims_h5, origin, spacing, 0_int32)
    end if
    deallocate(tmp)
  end subroutine save_r64_3d_wo_ghost



  ! ==========================
  ! === Core load helpers  ===
  ! ==========================
  subroutine load_core_int(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    integer(int32), allocatable, intent(out) :: A(..)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost

    integer(HID_T) :: file_id, dset_id, dspace_id
    integer        :: ierr, rank_h5
    integer(HSIZE_T), allocatable :: dims(:), maxdims(:)

    integer        :: myrank, nprocs, mpier
    integer        :: dims_i(3)      ! broadcast buffer (up to rank 3)
    integer        :: nelt, i
    logical        :: has_origin, has_spacing, has_nghost
    integer        :: ighost
    real(real64)   :: rtmp3(3)

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpier)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, mpier)

    ! Root reads from disk
    if (myrank == 0) then
      call h5open_f(ierr)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr); if (ierr/=0) error stop 'load_array_h5: h5fopen_f failed'
      call h5dopen_f(file_id, trim(dset), dset_id, ierr);              if (ierr/=0) error stop 'load_array_h5: h5dopen_f failed'
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sget_simple_extent_ndims_f(dspace_id, rank_h5, ierr)
      allocate(dims(rank_h5), maxdims(rank_h5))
      call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
    end if

    ! Broadcast rank and dims
    if (myrank == 0) then
      dims_i = 1
      do i=1, min(3, rank_h5)
        dims_i(i) = int(dims(i))
      end do
    end if
    call MPI_Bcast(rank_h5, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)
    call MPI_Bcast(dims_i, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)

    ! Allocate output array on non-root
    select rank (A)
    rank (1)
      if (myrank /= 0) allocate(A(dims_i(1)))
    rank (2)
      if (myrank /= 0) allocate(A(dims_i(1), dims_i(2)))
    rank (3)
      if (myrank /= 0) allocate(A(dims_i(1), dims_i(2), dims_i(3)))
    rank default
      error stop 'load_array_h5(int32): unexpected rank'
    end select

    ! Root reads dataset payload into A, then broadcast to all
    select rank (A)
    rank (1)
      if (myrank == 0) call h5dread_f(dset_id, H5T_NATIVE_INTEGER, A, dims, ierr)
      nelt = dims_i(1)
      call MPI_Bcast(A, nelt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)
    rank (2)
      if (myrank == 0) call h5dread_f(dset_id, H5T_NATIVE_INTEGER, A, dims, ierr)
      nelt = dims_i(1)*dims_i(2)
      call MPI_Bcast(A, nelt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)
    rank (3)
      if (myrank == 0) call h5dread_f(dset_id, H5T_NATIVE_INTEGER, A, dims, ierr)
      nelt = dims_i(1)*dims_i(2)*dims_i(3)
      call MPI_Bcast(A, nelt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)
    end select

    ! Broadcast optional attrs if the caller provided them (assume consistent presence on all ranks)
    has_origin  = present(origin)
    has_spacing = present(spacing)
    has_nghost  = present(nghost)
    call MPI_Bcast(has_origin,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpier)
    call MPI_Bcast(has_spacing, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpier)
    call MPI_Bcast(has_nghost,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpier)

    if (present(origin) .and. has_origin) then
      if (myrank == 0) then
        call read_attr_vec3(dset_id, "origin", rtmp3)
      end if
      call MPI_Bcast(rtmp3, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
      origin = rtmp3
    end if
    if (present(spacing) .and. has_spacing) then
      if (myrank == 0) then
        call read_attr_vec3(dset_id, "spacing", rtmp3)
      end if
      call MPI_Bcast(rtmp3, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
      spacing = rtmp3
    end if
    if (present(nghost) .and. has_nghost) then
      if (myrank == 0) then
        call read_attr_int1(dset_id, "nghost", ighost)
      end if
      call MPI_Bcast(ighost, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)
      nghost = ighost
    end if

    if (myrank == 0) then
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      deallocate(dims, maxdims)
    end if
  end subroutine load_core_int


  subroutine load_core_real(filename, dset, A, origin, spacing, nghost)
    character(*),               intent(in)  :: filename, dset
    real(real64),   allocatable, intent(out) :: A(..)
    real(real64),  dimension(3), intent(out), optional :: origin, spacing
    integer(int32),             intent(out), optional :: nghost

    integer(HID_T) :: file_id, dset_id, dspace_id
    integer        :: ierr, rank_h5
    integer(HSIZE_T), allocatable :: dims(:), maxdims(:)

    integer        :: myrank, nprocs, mpier
    integer        :: dims_i(3)
    integer        :: nelt, i
    logical        :: has_origin, has_spacing, has_nghost
    integer        :: ighost
    real(real64)   :: rtmp3(3)

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpier)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, mpier)

    if (myrank == 0) then
      call h5open_f(ierr)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr); if (ierr/=0) error stop 'load_array_h5: h5fopen_f failed'
      call h5dopen_f(file_id, trim(dset), dset_id, ierr);              if (ierr/=0) error stop 'load_array_h5: h5dopen_f failed'
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sget_simple_extent_ndims_f(dspace_id, rank_h5, ierr)
      allocate(dims(rank_h5), maxdims(rank_h5))
      call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
    end if

    if (myrank == 0) then
      dims_i = 1
      do i=1, min(3, rank_h5)
        dims_i(i) = int(dims(i))
      end do
    end if
    call MPI_Bcast(rank_h5, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)
    call MPI_Bcast(dims_i, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)

    select rank (A)
    rank (1)
      if (myrank /= 0) allocate(A(dims_i(1)))
    rank (2)
      if (myrank /= 0) allocate(A(dims_i(1), dims_i(2)))
    rank (3)
      if (myrank /= 0) allocate(A(dims_i(1), dims_i(2), dims_i(3)))
    rank default
      error stop 'load_array_h5(real64): unexpected rank'
    end select

    select rank (A)
    rank (1)
      if (myrank == 0) call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, ierr)
      nelt = dims_i(1)
      call MPI_Bcast(A, nelt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
    rank (2)
      if (myrank == 0) call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, ierr)
      nelt = dims_i(1)*dims_i(2)
      call MPI_Bcast(A, nelt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
    rank (3)
      if (myrank == 0) call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, ierr)
      nelt = dims_i(1)*dims_i(2)*dims_i(3)
      call MPI_Bcast(A, nelt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
    end select

    has_origin  = present(origin)
    has_spacing = present(spacing)
    has_nghost  = present(nghost)
    call MPI_Bcast(has_origin,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpier)
    call MPI_Bcast(has_spacing, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpier)
    call MPI_Bcast(has_nghost,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpier)

    if (present(origin) .and. has_origin) then
      if (myrank == 0) call read_attr_vec3(dset_id, "origin",  rtmp3)
      call MPI_Bcast(rtmp3, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
      origin = rtmp3
    end if
    if (present(spacing) .and. has_spacing) then
      if (myrank == 0) call read_attr_vec3(dset_id, "spacing", rtmp3)
      call MPI_Bcast(rtmp3, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
      spacing = rtmp3
    end if
    if (present(nghost) .and. has_nghost) then
      if (myrank == 0) call read_attr_int1(dset_id, "nghost", ighost)
      call MPI_Bcast(ighost, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpier)
      nghost = ighost
    end if

    if (myrank == 0) then
      call h5sclose_f(dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
      deallocate(dims, maxdims)
    end if
  end subroutine load_core_real


  ! --- small helpers (root-side) to read attributes safely ---
  subroutine read_attr_vec3(dset_id, name, outv)
    integer(HID_T), intent(in)  :: dset_id
    character(*),   intent(in)  :: name
    real(real64),   intent(out) :: outv(3)
    integer(HID_T) :: attr_id
    integer        :: ierr
    logical        :: has
    integer(HSIZE_T) :: three(1)
    three = [3_HSIZE_T]
    call h5aexists_f(dset_id, trim(name), has, ierr)
    if (has) then
      call h5aopen_f(dset_id, trim(name), attr_id, ierr)
      call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, outv, three, ierr)
      call h5aclose_f(attr_id, ierr)
    else
      outv = 0.0_real64
    end if
  end subroutine read_attr_vec3

  subroutine read_attr_int1(dset_id, name, outi)
    integer(HID_T), intent(in)  :: dset_id
    character(*),   intent(in)  :: name
    integer,        intent(out) :: outi
    integer(HID_T) :: attr_id
    integer        :: ierr
    logical        :: has
    integer(HSIZE_T) :: one(1)
    one = [1_HSIZE_T]
    call h5aexists_f(dset_id, trim(name), has, ierr)
    if (has) then
      call h5aopen_f(dset_id, trim(name), attr_id, ierr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, outi, one, ierr)
      call h5aclose_f(attr_id, ierr)
    else
      outi = 0
    end if
  end subroutine read_attr_int1


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
