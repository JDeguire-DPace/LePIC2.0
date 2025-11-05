module mod_grid
    use iso_fortran_env
    use mod_geometry
    implicit none
    private

    public :: Grid

    type :: Grid
        integer :: ncells_x1, ncells_x2, ncells_x3

        real(real64)  :: dx1, dx2, dx3, Length_x1, Length_x2, Length_x3
        integer :: nghost = 2
      contains
        procedure :: initialize_from_dom
        procedure :: total_cells_with_ghost_x1
        procedure :: total_cells_with_ghost_x2
        procedure :: total_cells_with_ghost_x3
    end type

    contains
    subroutine initialize_from_dom(self, dom)
        class(Grid), intent(inout) :: self
        type(Domain), intent(in) :: dom

        self%ncells_x1 = dom%ncell_x1; self%ncells_x2 = dom%ncell_x2; self%ncells_x3 = dom%ncell_x3
        self%Length_x1 = dom%Length_x1; self%Length_x2 = dom%Length_x2; self%Length_x3 = dom%Length_x3
        self%dx1 = dom%dx1; self%dx2 = dom%dx2; self%dx3 = dom%dx3

    end subroutine initialize_from_dom


    integer function total_cells_with_ghost_x1(self) result(nx1_total)
        class(Grid), intent(in) :: self

        nx1_total = self%ncells_x1 + 2*self%nghost
    end function total_cells_with_ghost_x1

    integer function total_cells_with_ghost_x2(self) result(nx2_total)
        class(Grid), intent(in) :: self
        nx2_total = self%ncells_x2 + 2*self%nghost
    end function total_cells_with_ghost_x2

    integer function total_cells_with_ghost_x3(self) result(nx3_total)
        class(Grid), intent(in) :: self
        nx3_total = self%ncells_x3 + 2*self%nghost
    end function total_cells_with_ghost_x3

end module mod_grid