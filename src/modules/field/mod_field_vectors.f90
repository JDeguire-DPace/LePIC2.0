! module mod_field_vectors
!   use iso_fortran_env, only: real64
!   use mod_grid       , only: Grid
!   use mod_parallel   , only: allreduce_sum_3d_real8
!   implicit none
!   private
!   public :: ElectricField, CurrentDensity, VectorPotential

!   type :: VectorField
!     real(real64), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
!   contains
!     procedure :: allocate_on
!     procedure :: zero
!   end type

!   type, extends(VectorField) :: ElectricField
!   contains
!     procedure :: compute_from_phi  ! central differences
!   end type

!   type, extends(VectorField) :: CurrentDensity
!   contains
!     procedure :: allreduce_sum => vf_allreduce
!   end type

!   type, extends(VectorField) :: VectorPotential
!   end type

! contains

!   subroutine allocate_on(f, g)
!     class(VectorField), intent(inout) :: f
!     type(Grid)        , intent(in)    :: g
!     allocate(f%x(0:g%nx+2,0:g%ny+2,0:g%nz+2))
!     allocate(f%y(0:g%nx+2,0:g%ny+2,0:g%nz+2))
!     allocate(f%z(0:g%nx+2,0:g%ny+2,0:g%nz+2))
!     call f%zero()
!   end subroutine

!   subroutine zero(f)
!     class(VectorField), intent(inout) :: f
!     f%x = 0.0_real64; f%y = 0.0_real64; f%z = 0.0_real64
!   end subroutine

!   subroutine compute_from_phi(E, phi, g)
!     class(ElectricField)     , intent(inout) :: E
!     type(ElectricPotential)  , intent(in)    :: phi
!     type(Grid)               , intent(in)    :: g
!     integer :: i,j,k
!     real(real64), parameter :: half=0.5_real64
! !$omp parallel do collapse(3)
!     do k=1,g%nz+1; do j=1,g%ny+1; do i=1,g%nx+1
!       E%x(i,j,k) = -(phi%a(i+1,j  ,k  ) - phi%a(i-1,j  ,k  ))*(half/g%dx)
!       E%y(i,j,k) = -(phi%a(i  ,j+1,k  ) - phi%a(i  ,j-1,k  ))*(half/g%dy)
!       E%z(i,j,k) = -(phi%a(i  ,j  ,k+1) - phi%a(i  ,j  ,k-1))*(half/g%dz)
!     end do; end do; end do
!   end subroutine

!   subroutine vf_allreduce(J)
!     class(CurrentDensity), intent(inout) :: J
!     call allreduce_sum_3d_real8(J%x)
!     call allreduce_sum_3d_real8(J%y)
!     call allreduce_sum_3d_real8(J%z)
!   end subroutine

! end module mod_field_vectors
