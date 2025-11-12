!==============================================================
!  mod_reactions.f90  (robust, ASCII-only I/O, no substring traps)
!  - Reads REACTION/REAC blocks & σ(E) tables from chemistry file
!  - σ(E) interpolation and (σ v_rel)_max ceilings
!  - Writes a summary into ./Outputs/particles.out
!  - Optional σ(E) dumps grouped by incident species
!==============================================================
module mod_reactions
  use iso_fortran_env, only: real64, int32, output_unit
  use iso_c_binding,   only: c_int
  use mod_InputFiles,  only: find_particles_file
  use mod_functionsText                 ! expects: is_dashed()
  use mod_particles,   only: SpeciesTable, species_index_of
  use mpi
  implicit none
  private

  !----------------------- Enumerated kinds -----------------------
  enum, bind(c)
    enumerator :: REACT_ELASTIC      = 1_c_int
    enumerator :: REACT_IONIZATION   = 2_c_int
    enumerator :: REACT_EXCITATION   = 3_c_int
    enumerator :: REACT_CHARGEEXCH   = 4_c_int
    enumerator :: REACT_DISSOCIATION = 5_c_int
    enumerator :: REACT_DATT         = 6_c_int
    enumerator :: REACT_DION         = 7_c_int
    enumerator :: REACT_DISS         = 8_c_int
    enumerator :: REACT_ELA          = 9_c_int
    enumerator :: REACT_EXC          = 10_c_int
    enumerator :: REACT_ION          = 11_c_int
    enumerator :: REACT_VIBB         = 12_c_int
    enumerator :: REACT_VIBS         = 13_c_int
    enumerator :: REACT_RRE          = 14_c_int
    enumerator :: REACT_DRE          = 15_c_int
    enumerator :: REACT_CEX          = 16_c_int
  end enum

  public :: REACT_ELASTIC, REACT_IONIZATION, REACT_EXCITATION, REACT_CHARGEEXCH, REACT_DISSOCIATION, &
            REACT_DATT, REACT_DION, REACT_DISS, REACT_ELA, REACT_EXC, REACT_ION, REACT_VIBB, REACT_VIBS, &
            REACT_RRE, REACT_DRE, REACT_CEX

  !---------------------- Public types / API ----------------------
  type :: CrossSectionTable
     real(real64), allocatable :: energy_eV_list(:)
     real(real64), allocatable :: cross_section_m2_list(:)
  end type

  type :: Reaction
     integer :: num_reactants = 0
     integer :: num_products  = 0
     integer :: reactant_species_index(4) = 0
     integer :: product_species_index(6)  = 0
     integer(c_int) :: type_id       = REACT_ELASTIC
     real(real64)   :: threshold_eV  = 0.0_real64
     real(real64)   :: heavy_gain_eV = 0.0_real64
     real(real64)   :: scale_energy  = 1.0_real64
     real(real64)   :: scale_sigma   = 1.0_real64
     integer        :: product_electron_count = 0
     type(CrossSectionTable) :: table
     character(len=:), allocatable :: label_text
  end type

  type :: ReactionSet
     type(Reaction), allocatable :: list(:)
  end type

  public :: CrossSectionTable, Reaction, ReactionSet
  public :: load_reactions_from_file
  public :: sigma_at_energy
  public :: count_reactions
  public :: compute_sigv_ceiling_all
  public :: write_reactions_to_particles_out
  public :: write_sigma_grouped_by_incident

contains

  !================================================================
  ! Read all REACTION/REAC blocks + σ(E) tables from chemistry file
  !================================================================
  subroutine load_reactions_from_file(out_reactions, species_catalog)
    type(ReactionSet), intent(out) :: out_reactions
    type(SpeciesTable),intent(in)  :: species_catalog
    integer :: u, ios

    if (allocated(out_reactions%list)) deallocate(out_reactions%list)
    allocate(out_reactions%list(0))

    open(newunit=u, file=trim(find_particles_file()), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'mod_reactions: cannot open chemistry file.'

    call parse_all_reactions(u, out_reactions, species_catalog)
    close(u)

    write(output_unit,'(A,I0)') 'Reactions loaded = ', count_reactions(out_reactions)
    call flush(output_unit)
  end subroutine load_reactions_from_file

  subroutine parse_all_reactions(u, out_reactions, species_catalog)
    integer,             intent(in)    :: u
    type(ReactionSet),   intent(inout) :: out_reactions
    type(SpeciesTable),  intent(in)    :: species_catalog
    character(len=512) :: line
    character(len=:), allocatable :: t
    integer :: ios

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) exit
      call strip_comment(line)
      t = adjustl(trim(line))
      if (len_trim(t) == 0) cycle

      ! Avoid char-array literal to dodge substring bugs on some compilers
      if (starts_with_ci(t,'REACTION') .or. starts_with_ci(t,'REAC')) then
        call parse_one_reaction_block(u, out_reactions, species_catalog)
      end if
    end do
  end subroutine parse_all_reactions

  subroutine parse_one_reaction_block(u, out_reactions, species_catalog)
    integer,             intent(in)    :: u
    type(ReactionSet),   intent(inout) :: out_reactions
    type(SpeciesTable),  intent(in)    :: species_catalog

    type(Reaction) :: rx
    integer :: ios
    character(len=512) :: line

    ! 1) Label
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing reaction label.')
    call strip_comment(line); rx%label_text = trim(line)
    call parse_label_species(rx%label_text, rx, species_catalog)

    ! 2) Eth  dE
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing "Eth dE".'); call strip_comment(line)
    call read_two_reals(line, rx%threshold_eV, rx%heavy_gain_eV)

    ! 3) scale_E  scale_sigma
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing "scale_E scale_sigma".'); call strip_comment(line)
    call read_two_reals(line, rx%scale_energy, rx%scale_sigma)

    ! 4) Kind
    read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Missing reaction kind.'); call strip_comment(line)
    rx%type_id = map_kind_string(trim(line))

    ! 5) Seek σ(E) table start
    do
      read(u,'(A)', iostat=ios) line; call expect_ok(ios, 'Unexpected EOF before sigma(E) table.')
      call strip_comment(line)
      if (len_trim(line) == 0) cycle
      if (is_dashed(trim(line))) then
        exit
      else
        block
          real(real64) :: etry, stry
          integer :: ios2
          read(line,*,iostat=ios2) etry, stry
          if (ios2 == 0) then
            backspace(u)
            exit
          end if
        end block
      end if
    end do

    ! 6) σ(E) pairs to next dashed line
    call read_sigma_pairs(u, rx%table, rx%scale_energy, rx%scale_sigma)

    ! 7) Count product electrons
    rx%product_electron_count = count(rx%product_species_index(1:rx%num_products) == &
                                      species_index_of(species_catalog, '[e]'))

    ! 8) Store
    call push_reaction(out_reactions, rx)
  end subroutine parse_one_reaction_block

  !================================================================
  ! σ(E) interpolation (linear, clamped)
  !================================================================
  pure function sigma_at_energy(table, energy_eV) result(sigma_m2)
    type(CrossSectionTable), intent(in) :: table
    real(real64),            intent(in) :: energy_eV
    real(real64) :: sigma_m2
    integer :: n, iL

    sigma_m2 = 0.0_real64
    n = size(table%energy_eV_list); if (n == 0) return

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
  ! Ceilings for null-collision method
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

  pure function sigv_ceiling_one(rx, species_catalog) result(max_sigv)
    type(Reaction),     intent(in) :: rx
    type(SpeciesTable), intent(in) :: species_catalog
    real(real64) :: max_sigv, mu, E, vr, m1, m2
    integer :: n, i
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

  !================================================================
  ! σ(E) table reader
  !================================================================
  subroutine read_sigma_pairs(u, table, scale_E, scale_sigma)
    integer,                 intent(in)  :: u
    type(CrossSectionTable), intent(out) :: table
    real(real64),            intent(in)  :: scale_E, scale_sigma

    character(len=256) :: line
    integer :: ios
    real(real64) :: e_val, s_val

    if (allocated(table%energy_eV_list)) then
      deallocate(table%energy_eV_list, table%cross_section_m2_list)
    end if
    allocate(table%energy_eV_list(0), table%cross_section_m2_list(0))

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) exit
      call strip_comment(line)
      if (len_trim(line) == 0) cycle
      if (is_dashed(trim(line))) exit

      if (.not. looks_like_two_reals(line)) cycle
      read(line,*,iostat=ios) e_val, s_val
      if (ios == 0) call push_sigma_point(table, e_val*scale_E, max(0.0_real64, s_val*scale_sigma))
    end do

    call enforce_strictly_increasing(table%energy_eV_list)
  end subroutine read_sigma_pairs

  !================================================================
  ! Species parsing in labels like: "[e]+[D2]->[D-]+[D]"
  !================================================================
  subroutine parse_label_species(label, rx, species_catalog)
    character(len=*), intent(in)    :: label
    type(Reaction),   intent(inout) :: rx
    type(SpeciesTable),intent(in)   :: species_catalog

    character(len=:), allocatable :: lhs, rhs
    character(len=256) :: token_list(128)
    integer :: n_tok, i, idx, ierr, rank0

    call MPI_Comm_rank(MPI_COMM_WORLD, rank0, ierr)

    call split_by_arrow(trim(label), lhs, rhs)

    call split_species_list(lhs, token_list, n_tok)
    rx%num_reactants = 0
    do i=1,n_tok
      idx = species_index_of(species_catalog, trim(token_list(i)))
      if (idx > 0) then
        if (rx%num_reactants < size(rx%reactant_species_index)) then
          rx%num_reactants = rx%num_reactants + 1
          rx%reactant_species_index(rx%num_reactants) = idx
        else
          if (rank0 == 0) write(output_unit,'(A)') 'WARNING: too many reactants (extra ignored): '//trim(label)
        end if
      else
        if (rank0 == 0) write(output_unit,'(A)') 'Unknown reactant: '//trim(token_list(i))//'  <<'//trim(label)//'>>'
      end if
    end do

    call split_species_list(rhs, token_list, n_tok)
    rx%num_products = 0
    do i=1,n_tok
      idx = species_index_of(species_catalog, trim(token_list(i)))
      if (idx > 0) then
        if (rx%num_products < size(rx%product_species_index)) then
          rx%num_products = rx%num_products + 1
          rx%product_species_index(rx%num_products) = idx
        else
          if (rank0 == 0) write(output_unit,'(A)') 'WARNING: too many products (extra ignored): '//trim(label)
        end if
      else
        if (rank0 == 0) write(output_unit,'(A)') 'Unknown product: '//trim(token_list(i))//'  <<'//trim(label)//'>>'
      end if
    end do
  end subroutine parse_label_species

  !---------------------- small helpers -------------------------
  pure logical function looks_like_two_reals(s) result(ok)
    character(len=*), intent(in) :: s
    real(real64) :: a, b
    integer :: ios
    ok = .false.
    read(s, *, iostat=ios) a, b
    if (ios == 0) ok = .true.
  end function looks_like_two_reals

  pure function lin_interp(x, xL, xR, yL, yR) result(y)
    real(real64), intent(in) :: x, xL, xR, yL, yR
    real(real64) :: y
    if (xR > xL) then
      y = yL + (x - xL) * (yR - yL) / (xR - xL)
    else
      y = yL
    end if
  end function lin_interp

  subroutine push_sigma_point(table, e_val, s_val)
    type(CrossSectionTable), intent(inout) :: table
    real(real64),            intent(in)    :: e_val, s_val
    integer :: n
    real(real64), allocatable :: ebuf(:), cbuf(:)
    n = size(table%energy_eV_list)
    allocate(ebuf(n+1), cbuf(n+1))
    if (n > 0) then
      ebuf(1:n) = table%energy_eV_list
      cbuf(1:n) = table%cross_section_m2_list
    end if
    ebuf(n+1) = e_val
    cbuf(n+1) = s_val
    if (allocated(table%energy_eV_list)) deallocate(table%energy_eV_list)
    if (allocated(table%cross_section_m2_list)) deallocate(table%cross_section_m2_list)
    call move_alloc(ebuf, table%energy_eV_list)
    call move_alloc(cbuf, table%cross_section_m2_list)
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
    case ("ELASTIC")     ; k = REACT_ELASTIC
    case ("IONIZATION")  ; k = REACT_IONIZATION
    case ("EXCITATION")  ; k = REACT_EXCITATION
    case ("CHARGEEXCHANGE"); k = REACT_CHARGEEXCH
    case ("DISSOCIATION"); k = REACT_DISSOCIATION
    case ("DATT")        ; k = REACT_DATT
    case ("DION")        ; k = REACT_DION
    case ("DISS")        ; k = REACT_DISS
    case ("ELA")         ; k = REACT_ELA
    case ("EXC")         ; k = REACT_EXC
    case ("ION")         ; k = REACT_ION
    case ("VIBB")        ; k = REACT_VIBB
    case ("VIBS")        ; k = REACT_VIBS
    case ("RRE")         ; k = REACT_RRE
    case ("DRE")         ; k = REACT_DRE
    case ("CEX")         ; k = REACT_CEX
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

  pure logical function starts_with_ci(s, prefix)
    character(len=*), intent(in) :: s, prefix
    integer :: ns, np
    ns = len_trim(s); np = len_trim(prefix)
    if (ns < np) then
      starts_with_ci = .false.
    else
      starts_with_ci = (uppercase(s(:np)) == uppercase(prefix(:np)))
    end if
  end function starts_with_ci

  subroutine split_species_list(side, tokens, n_tok)
    character(len=*), intent(in)  :: side
    character(len=*), intent(out) :: tokens(:)
    integer,          intent(out) :: n_tok
    character(len=:), allocatable :: s
    integer :: i, n, start_pos, depth, cap

    s = adjustl(trim(side))
    n = len_trim(s)
    n_tok = 0
    start_pos = 1
    depth = 0
    cap = size(tokens)

    do i = 1, max(1,n)
       select case (s(i:i))
       case ('['); depth = depth + 1
       case (']'); if (depth > 0) depth = depth - 1
       case ('+')
          if (depth == 0) then
             if (i-1 >= start_pos) then
                if (n_tok < cap) then
                  n_tok = n_tok + 1
                  tokens(n_tok) = adjustl(trim(s(start_pos:i-1)))
                end if
             end if
             start_pos = i + 1
          end if
       end select
    end do

    if (n >= start_pos) then
       if (n_tok < cap) then
         n_tok = n_tok + 1
         tokens(n_tok) = adjustl(trim(s(start_pos:n)))
       end if
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

  !==============================================================
  !  Writers
  !==============================================================
  subroutine write_reactions_to_particles_out(reactions, species_catalog, write_sigma_files)
    type(ReactionSet), intent(in) :: reactions
    type(SpeciesTable),intent(in) :: species_catalog
    logical, optional, intent(in) :: write_sigma_files

    integer :: ierr, rank_local, i, u, ios
    character(len=*), parameter :: outFile = './Outputs/particles.out'
    logical :: do_sigma
    do_sigma = .false.; if (present(write_sigma_files)) do_sigma = write_sigma_files

    call MPI_Comm_rank(MPI_COMM_WORLD, rank_local, ierr)
    if (rank_local /= 0) return

    call execute_command_line('mkdir -p ./Outputs', exitstat=ios)
    open(newunit=u, file=outFile, status='unknown', action='write', position='append', iostat=ios)
    if (ios /= 0) then
      write(output_unit,'(A)') 'write_reactions_to_particles_out: cannot open ./Outputs/particles.out'
      call flush(output_unit)
      return
    end if

    write(u,'(A)') ' '
    write(u,'(A)') 'Reactions:'
    do i = 1, count_reactions(reactions)
      write(u,'(I6,1X,A,1X,"        Eth=",ES12.4,1X,"        type=",A)') i, &
        trim(default_if_empty(reactions%list(i)%label_text,'(unlabeled)')), &
        reactions%list(i)%threshold_eV, trim(reaction_type_name(reactions%list(i)%type_id))
    end do
    write(u,'(A,I0)') 'Total number of reactions: ', count_reactions(reactions)
    close(u)

    if (do_sigma) call write_sigma_grouped_by_incident(reactions, species_catalog)
  end subroutine write_reactions_to_particles_out

  subroutine write_sigma_grouped_by_incident(reactions, species_catalog)
    type(ReactionSet), intent(in) :: reactions
    type(SpeciesTable),intent(in) :: species_catalog
    integer :: ierr, rank_local, i, s, u, ios, n, k
    character(len=32)  :: idxstr
    character(len=256) :: path

    call MPI_Comm_rank(MPI_COMM_WORLD, rank_local, ierr)
    if (rank_local /= 0) return

    call execute_command_line('mkdir -p ./Outputs', exitstat=ios)

    do s = 1, species_catalog%number_species
      write(idxstr,'(I0)') s
      path = './Outputs/reactions.'//trim(idxstr)
      open(newunit=u, file=path, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
        write(output_unit,'(A)') 'write_sigma_grouped_by_incident: cannot open file'
        call flush(output_unit)
        cycle
      end if

      do i = 1, count_reactions(reactions)
        if (reactions%list(i)%num_reactants >= 1) then
          if (reactions%list(i)%reactant_species_index(1) == s) then
            n = size(reactions%list(i)%table%energy_eV_list)
            if (n > 0) then
              write(u,'(A)') trim(default_if_empty(reactions%list(i)%label_text,'(unlabeled)'))
              write(u,'(A)') '# E[eV]    sigma[m^2]'
              do k = 1, n
                write(u,'(F12.4,1X,ES14.6)') reactions%list(i)%table%energy_eV_list(k), &
                                            reactions%list(i)%table%cross_section_m2_list(k)
              end do
              write(u,'(A)') ' '
            end if
          end if
        end if
      end do
      close(u)
    end do
  end subroutine write_sigma_grouped_by_incident

  pure function reaction_type_name(k) result(name)
    integer(c_int), intent(in) :: k
    character(len=24) :: name
    select case (k)
    case (REACT_ELASTIC)     ; name = 'ELASTIC'
    case (REACT_IONIZATION)  ; name = 'IONIZATION'
    case (REACT_EXCITATION)  ; name = 'EXCITATION'
    case (REACT_CHARGEEXCH)  ; name = 'CHARGEEXCHANGE'
    case (REACT_DISSOCIATION); name = 'DISSOCIATION'
    case (REACT_DATT)        ; name = 'DATT'
    case (REACT_DION)        ; name = 'DION'
    case (REACT_DISS)        ; name = 'DISS'
    case (REACT_ELA)         ; name = 'ELA'
    case (REACT_EXC)         ; name = 'EXC'
    case (REACT_ION)         ; name = 'ION'
    case (REACT_VIBB)        ; name = 'VIBB'
    case (REACT_VIBS)        ; name = 'VIBS'
    case (REACT_RRE)         ; name = 'RRE'
    case (REACT_DRE)         ; name = 'DRE'
    case (REACT_CEX)         ; name = 'CEX'
    case default             ; name = 'UNKNOWN'
    end select
  end function reaction_type_name

  pure function default_if_empty(s, fallback) result(out)
    character(len=*), intent(in) :: s, fallback
    character(len=:), allocatable :: out
    if (len_trim(s) == 0) then
      out = trim(fallback)
    else
      out = trim(s)
    end if
  end function default_if_empty

  pure integer function count_reactions(reactions) result(n)
    type(ReactionSet), intent(in) :: reactions
    if (.not. allocated(reactions%list)) then
      n = 0
    else
      n = size(reactions%list)
    end if
  end function count_reactions

end module mod_reactions