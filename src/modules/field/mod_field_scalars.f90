! module mod_field_scalars
!   use iso_fortran_env, only: real64
!   use mod_grid       , only: Grid
!   use mod_parallel   , only: par_is_root, allreduce_sum_3d_real8, bcast_3d_real8
!   implicit none
!   private
!   public :: ElectricPotential, ChargeDensity
!   public :: ScalarField       ! base you can reuse for Ax,Ay,Az, etc.

!   type :: ScalarField
!     real(real64), allocatable :: a(:,:,:)
!   contains
!     procedure :: allocate_on
!     procedure :: zero
!     procedure :: data  => sf_data
!   end type

!   type, extends(ScalarField) :: ElectricPotential
!   contains
!     procedure :: bcast => sf_bcast
!   end type

!   type, extends(ScalarField) :: ChargeDensity
!   contains
!     procedure :: allreduce_sum => sf_allreduce
!   end type

!   contains

!   subroutine allocate_on(f, g)
!     class(ScalarField), intent(inout) :: f
!     type(Grid)        , intent(in)    :: g
!     if (allocated(f%a)) deallocate(f%a)
!     allocate(f%a(0:g%nx+2, 0:g%ny+2, 0:g%nz+2))
!     f%a = 0.0_real64
!   end subroutine

!   subroutine zero(f)
!     class(ScalarField), intent(inout) :: f
!     f%a = 0.0_real64
!   end subroutine

!   function sf_data(f) result(p)
!     class(ScalarField), target, intent(inout) :: f
!     real(real64), pointer :: p(:,:,:)
!     p => f%a
!   end function

!   subroutine sf_allreduce(f)
!     class(ChargeDensity), intent(inout) :: f
!     call allreduce_sum_3d_real8(f%a)
!   end subroutine

!   subroutine sf_bcast(f, root)
!     class(ElectricPotential), intent(inout) :: f
!     integer, intent(in) :: root
!     call bcast_3d_real8(f%a, root)
!   end subroutine

! end module mod_field_scalars
