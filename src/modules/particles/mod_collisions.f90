!==============================================================
!  mod_collisions.f90  (clean, OOP-friendly, insightful names)
!  - Reads REACTION blocks & σ(E) tables from chemistry file
!  - Provides σ(E) interpolation and (σ v_rel)_max ceilings
!  - Uses proper Fortran ENUM with iso_c_binding (c_int)
!==============================================================
module mod_collisions
  use iso_fortran_env, only: real64, int32, output_unit
  use iso_c_binding,   only: c_int
  use mod_InputFiles,  only: find_particles_file
  use mod_functionsText            ! expecting is_dashed(), etc.
  use mod_particles,   only: SpeciesTable, species_index_of
  implicit none
  private

  !----------------------- Enumerated kinds -----------------------
  ! NOTE: export *enumerators*; there is no named enum type.
  enum, bind(c)
    enumerator :: REACT_ELASTIC      = 1_c_int
    enumerator :: REACT_IONIZATION   = 2_c_int
    enumerator :: REACT_EXCITATION   = 3_c_int
    enumerator :: REACT_CHARGEEXCH   = 4_c_int
    enumerator :: REACT_DISSOCIATION = 5_c_int
    ! extra short-tags you mentioned
    enumerator :: REACT_DATT         = 6_c_int
    enumerator :: REACT_DION         = 7_c_int
    enumerator :: REACT_DISS         = 8_c_int
    enumerator :: REACT_ELA          = 9_c_int
    enumerator :: REACT_EXC          = 10_c_int
  end enum

  public :: REACT_ELASTIC, REACT_IONIZATION, REACT_EXCITATION, REACT_CHARGEEXCH, REACT_DISSOCIATION, &
            REACT_DATT, REACT_DION, REACT_DISS, REACT_ELA, REACT_EXC

  !---------------------- Public types / API ----------------------
  public :: CrossSectionTable, Reaction, ReactionSet
  public :: load_reactions_from_file
  public :: sigma_at_energy
  public :: count_reactions, reaction_at
  public :: compute_sigv_ceiling_all

  ! Tabulated σ(E) (energies strictly increasing; units: eV, m²)
  type :: CrossSectionTable
     real(real64), allocatable :: energy_eV_list(:)
     real(real64), allocatable :: cross_section_m2_list(:)
  end type

  ! One reaction with species indices resolved against SpeciesTable
  type :: Reaction
     ! Species mapping
     integer :: num_reactants = 0
     integer :: num_products  = 0
     integer :: reactant_species_index(4) = 0  !! up to 4 reactants
     integer :: product_species_index(6)  = 0  !! up to 6 products

     ! Kind and metadata
     integer(c_int) :: type_id       = REACT_ELASTIC
     real(real64)   :: threshold_eV  = 0.0_real64   !! Eth
     real(real64)   :: heavy_gain_eV = 0.0_real64   !! dE to heavy byproducts
     real(real64)   :: scale_energy  = 1.0_real64   !! scale applied to E column
     real(real64)   :: scale_sigma   = 1.0_real64   !! scale applied to σ column
     integer        :: product_electron_count = 0

     type(CrossSectionTable) :: table
     character(len=:), allocatable :: label_text   !! "[e]+[H]->[e]+[H+]+[e]"
  end type

  type :: ReactionSet
     type(Reaction), allocatable :: list(:)
  end type

  contains
  !================================================================
  ! Read all REACTION blocks + σ(E) tables from chemistry file
  !================================================================
  subroutine load_reactions_from_file(out_reactions, species_catalog)
    type(ReactionSet), intent(out) :: out_reactions
    type(SpeciesTable),intent(in)  :: species_catalog
    integer :: u, ios
    character(len=512) :: line

    if (allocated(out_reactions%list)) deallocate(out_reactions%list)
    allocate(out_reactions%list(0))

    open(newunit=u, file=trim(find_particles_file()), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'mod_collisions: cannot open chemistry file.'

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) exit
      call strip_comment(line)
      if (len_trim(line) == 0) cycle

      if (starts_with(trim(line), 'REACTION')) then
        call parse_one_reaction_block(u, out_reactions, species_catalog)
      end if
    end do
    close(u)

    write(output_unit,'(A,I0)') 'Collisions: reactions loaded = ', count_reactions(out_reactions)
  end subroutine load_reactions_from_file

  !================================================================
  ! σ(E) interpolation (linear, clamped at bounds)
  !================================================================
  pure function sigma_at_energy(table, energy_eV) result(sigma_m2)
    type(CrossSectionTable), intent(in) :: table
    real(real64),            intent(in) :: energy_eV
    real(real64) :: sigma_m2
    integer :: n, iL

    sigma_m2 = 0.0_real64
    n = size(table%energy_eV_list)
    if (n == 0) return

    if (energy_eV <= table%energy_eV_list(1)) then
      sigma_m2 = table%cross_section_m2_list(1); return
    end if
    if (energy_eV >= table%energy_eV_list(n)) then
      sigma_m2 = table%cross_section_m2_list(n); return
    end if

    do iL = 1, n-1
      if (energy_eV >= table%energy_eV_list(iL) .and. energy_eV <= table%energy_eV_list(iL+1)) then
        sigma_m2 = lin_interp(energy_eV,                             &
                              table%energy_eV_list(iL),               &
                              table%energy_eV_list(iL+1),             &
                              table%cross_section_m2_list(iL),        &
                              table%cross_section_m2_list(iL+1))
        return
      end if
    end do
  end function sigma_at_energy

  !================================================================
  ! Utilities for sets
  !================================================================
  pure function count_reactions(reactions) result(n)
    type(ReactionSet), intent(in) :: reactions
    integer :: n
    if (.not. allocated(reactions%list)) then
      n = 0
    else
      n = size(reactions%list)
    end if
  end function count_reactions

  pure function reaction_at(reactions, i) result(r)
    type(ReactionSet), intent(in) :: reactions
    integer,           intent(in) :: i
    type(Reaction) :: r
    r = reactions%list(i)
  end function reaction_at

  !================================================================
  ! (σ v_rel)_max per reaction (null-collision ceiling)
  !   v_rel = sqrt(2 * E * e / μ), μ = reduced mass of first 2 reactants
  !================================================================
  subroutine compute_sigv_ceiling_all(reactions, species_catalog, sigv_ceiling)
    type(ReactionSet), intent(in)  :: reactions
    type(SpeciesTable),intent(in)  :: species_catalog
    real(real64), allocatable, intent(out) :: sigv_ceiling(:)
    integer :: i, n

    n = count_reactions(reactions)
    allocate(sigv_ceiling(n)); sigv_ceiling = 0.0_real64

    do i = 1, n
      sigv_ceiling(i) = sigv_ceiling_one(reactions%list(i), species_catalog)
    end do
  end subroutine compute_sigv_ceiling_all

  !======================== Internal ============================

  subroutine parse_one_reaction_block(u, out_reactions, species_catalog)
    integer,             intent(in)    :: u
    type(ReactionSet),   intent(inout) :: out_reactions
    type(SpeciesTable),  intent(in)    :: species_catalog

    type(Reaction) :: rx
    integer :: ios
    character(len=512) :: line

    ! label
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing reaction label.')
    call strip_comment(line); rx%label_text = trim(line)
    call parse_label_species(rx%label_text, rx, species_catalog)

    ! Eth dE
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing "Eth dE".'); call strip_comment(line)
    call read_two_reals(line, rx%threshold_eV, rx%heavy_gain_eV)

    ! scale_E scale_sigma
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing "scale_E scale_sigma".'); call strip_comment(line)
    call read_two_reals(line, rx%scale_energy, rx%scale_sigma)

    ! kind
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing reaction kind.'); call strip_comment(line)
    rx%type_id = map_kind_string(trim(line))

    ! seek NOTES
    do
      read(u,'(A)', iostat=ios) line; if (ios /= 0) exit
      call strip_comment(line)
      if (index(uppercase(trim(line)), 'NOTES') > 0) exit
    end do
    ! skip header line
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Unexpected EOF before σ(E) table.')

    ! read E-σ pairs
    call read_sigma_pairs(u, rx%table, rx%scale_energy, rx%scale_sigma)

    ! count product electrons
    rx%product_electron_count = count(rx%product_species_index(1:rx%num_products) == &
                                      species_index_of(species_catalog, '[e]'))

    call push_reaction(out_reactions, rx)
  end subroutine parse_one_reaction_block

  pure function sigv_ceiling_one(rx, species_catalog) result(max_sigv)
    type(Reaction),     intent(in) :: rx
    type(SpeciesTable), intent(in) :: species_catalog
    real(real64) :: max_sigv
    real(real64) :: mu, E, vr
    integer :: n, i
    real(real64) :: m1, m2
    real(real64), parameter :: e_charge = 1.602176634e-19_real64

    max_sigv = 0.0_real64
    if (rx%num_reactants < 2) return

    m1 = abs(species_catalog%species(rx%reactant_species_index(1))%mass)
    m2 = abs(species_catalog%species(rx%reactant_species_index(2))%mass)
    if (m1 <= 0.0_real64 .or. m2 <= 0.0_real64) return

    mu = (m1*m2)/(m1+m2)
    n  = size(rx%table%energy_eV_list)

    do i = 1, n
      E  = max(1.0_real64, rx%table%energy_eV_list(i))  ! avoid 0 eV
      vr = sqrt(2.0_real64 * E * e_charge / mu)
      max_sigv = max(max_sigv, vr * rx%table%cross_section_m2_list(i))
    end do
  end function sigv_ceiling_one

  subroutine read_sigma_pairs(u, table, scale_E, scale_sigma)
    integer,                intent(in)  :: u
    type(CrossSectionTable),intent(out) :: table
    real(real64),           intent(in)  :: scale_E, scale_sigma
    character(len=256) :: line
    integer :: ios
    real(real64) :: E, S

    if (allocated(table%energy_eV_list)) deallocate(table%energy_eV_list, table%cross_section_m2_list)
    allocate(table%energy_eV_list(0), table%cross_section_m2_list(0))

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) exit
      call strip_comment(line)
      if (len_trim(line) == 0) cycle
      if (is_dashed(trim(line))) exit

      read(line,*,iostat=ios) E, S
      if (ios /= 0) cycle
      call push_sigma_point(table, E*scale_E, max(0.0_real64, S*scale_sigma))
    end do

    call enforce_strictly_increasing(table%energy_eV_list)
  end subroutine read_sigma_pairs

  subroutine parse_label_species(label, rx, species_catalog)
    character(len=*), intent(in)    :: label
    type(Reaction),   intent(inout) :: rx
    type(SpeciesTable),intent(in)   :: species_catalog

    character(len=:), allocatable :: lhs, rhs
    character(len=256) :: token_list(32)
    integer :: n_tok, i, idx

    call split_by_arrow(trim(label), lhs, rhs)

    ! reactants (split '+' only outside brackets)
    call split_species_list(lhs, token_list, n_tok)
    rx%num_reactants = 0
    do i=1,n_tok
      idx = species_index_of(species_catalog, trim(token_list(i)))
      if (idx > 0) then
        rx%num_reactants = rx%num_reactants + 1
        rx%reactant_species_index(rx%num_reactants) = idx
      else
        write(output_unit,'(A)') 'Unknown reactant: '//trim(token_list(i))//'  <<'//trim(label)//'>>'
      end if
    end do

    ! products (split '+' only outside brackets)
    call split_species_list(rhs, token_list, n_tok)
    rx%num_products = 0
    do i=1,n_tok
      idx = species_index_of(species_catalog, trim(token_list(i)))
      if (idx > 0) then
        rx%num_products = rx%num_products + 1
        rx%product_species_index(rx%num_products) = idx
      else
        write(output_unit,'(A)') 'Unknown product: '//trim(token_list(i))//'  <<'//trim(label)//'>>'
      end if
    end do
  end subroutine parse_label_species

  !---------------------- small helpers -------------------------
  pure function lin_interp(x, xL, xR, yL, yR) result(y)
    real(real64), intent(in) :: x, xL, xR, yL, yR
    real(real64) :: y
    if (xR > xL) then
      y = yL + (x - xL) * (yR - yL) / (xR - xL)
    else
      y = yL
    end if
  end function lin_interp

  subroutine push_sigma_point(table, E, S)
    type(CrossSectionTable), intent(inout) :: table
    real(real64),            intent(in)    :: E, S
    integer :: n
    real(real64), allocatable :: e_tmp(:), s_tmp(:)
    n = size(table%energy_eV_list)
    allocate(e_tmp(n+1), s_tmp(n+1))
    if (n > 0) then
      e_tmp(1:n) = table%energy_eV_list
      s_tmp(1:n) = table%cross_section_m2_list
      deallocate(table%energy_eV_list, table%cross_section_m2_list)
    end if
    e_tmp(n+1) = E
    s_tmp(n+1) = S
    call move_alloc(e_tmp, table%energy_eV_list)
    call move_alloc(s_tmp, table%cross_section_m2_list)
  end subroutine push_sigma_point

  subroutine enforce_strictly_increasing(a)
    real(real64), intent(inout) :: a(:)
    integer :: i, n
    real(real64) :: last, eps
    n = size(a); if (n <= 1) return
    eps = epsilon(1.0_real64)
    last = a(1)
    do i = 2, n
      if (a(i) <= last) a(i) = last + max(eps, abs(last)*eps)
      last = a(i)
    end do
  end subroutine enforce_strictly_increasing

  integer(c_int) function map_kind_string(s) result(k)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: u
    u = uppercase(trim(s))
    select case (u)
    case ('COLLISION');      k = REACT_ELASTIC
    case ('IONIZATION');     k = REACT_IONIZATION
    case ('EXCITATION');     k = REACT_EXCITATION
    case ('CHARGEEXCHANGE'); k = REACT_CHARGEEXCH
    case ('DISSOCIATION');   k = REACT_DISSOCIATION
    case ('DATT');           k = REACT_DATT
    case ('DION');           k = REACT_DION
    case ('DISS');           k = REACT_DISS
    case ('ELA');            k = REACT_ELA
    case ('EXC');            k = REACT_EXC
    case default;            k = REACT_ELASTIC
    end select
  end function map_kind_string

  pure function uppercase(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len_trim(s)) :: t
    integer :: i, n, a
    n = len_trim(s)
    do i=1,n
      a = iachar(s(i:i))
      if (a >= iachar('a') .and. a <= iachar('z')) then
        t(i:i) = achar(a - 32)
      else
        t(i:i) = s(i:i)
      end if
    end do
  end function uppercase

  ! Split on '+' only when not inside [...]
  subroutine split_species_list(side, tokens, n_tok)
    character(len=*), intent(in)  :: side
    character(len=*), intent(out) :: tokens(:)
    integer,          intent(out) :: n_tok

    character(len=:), allocatable :: s
    integer :: i, n, start_pos, depth
    character(len=1) :: ch

    s = adjustl(trim(side))
    n = len_trim(s)
    n_tok = 0
    start_pos = 1
    depth = 0  ! bracket nesting depth

    do i = 1, max(1,n)
       ch = s(i:i)
       select case (ch)
       case ('['); depth = depth + 1
       case (']'); if (depth > 0) depth = depth - 1
       case ('+')
          if (depth == 0) then
             if (i-1 >= start_pos) then
                n_tok = n_tok + 1
                tokens(n_tok) = adjustl(trim(s(start_pos:i-1)))
             end if
             start_pos = i + 1
          end if
       end select
    end do

    if (n >= start_pos) then
       n_tok = n_tok + 1
       tokens(n_tok) = adjustl(trim(s(start_pos:n)))
    end if
  end subroutine split_species_list

  subroutine split_by_arrow(label, lhs, rhs)
    character(len=*), intent(in)  :: label
    character(len=:), allocatable, intent(out) :: lhs, rhs
    integer :: pos
    pos = index(label, '->')
    if (pos <= 0) then
      lhs = trim(label); rhs = ''
    else
      lhs = trim(label(:pos-1))
      rhs = trim(label(pos+2:))
    end if
  end subroutine split_by_arrow

  pure logical function starts_with(s, prefix)
    character(len=*), intent(in) :: s, prefix
    integer :: n
    n = len_trim(prefix)
    if (len_trim(s) < n) then
      starts_with = .false.
    else
      starts_with = (s(1:n) == prefix(1:n))
    end if
  end function starts_with

  subroutine strip_comment(line)
    character(len=*), intent(inout) :: line
    integer :: p
    p = index(line, '!')
    if (p > 0) line = line(:p-1)
  end subroutine strip_comment

  subroutine read_two_reals(text, a, b)
    character(len=*), intent(in)  :: text
    real(real64),     intent(out) :: a, b
    integer :: ios
    read(text,*,iostat=ios) a, b
    call expect_ok(ios, 'Expected two real numbers.')
  end subroutine read_two_reals

  subroutine push_reaction(set, rx)
    type(ReactionSet), intent(inout) :: set
    type(Reaction),    intent(in)    :: rx
    type(Reaction), allocatable :: tmp(:)
    integer :: n
    n = count_reactions(set)
    allocate(tmp(n+1))
    if (n > 0) tmp(1:n) = set%list
    tmp(n+1) = rx
    if (allocated(set%list)) deallocate(set%list)
    call move_alloc(tmp, set%list)
  end subroutine push_reaction

  subroutine expect_ok(ios, msg)
    integer,         intent(in) :: ios
    character(len=*),intent(in) :: msg
    if (ios /= 0) error stop trim(msg)
  end subroutine expect_ok

end module mod_collisions
