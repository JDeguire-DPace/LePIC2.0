!==============================================================
!  mod_particles.f90
!==============================================================
module mod_particles
  !! Particle containers + species metadata + builders
  use iso_fortran_env, only: real64, int32, int8, int64
  use mpi
  use mod_InputFiles
  use mod_functionsText          ! expects split_ws, etc.
  use mod_constants              ! expects: elementary_charge, amu, mass_electron
  use mod_geometry
  implicit none
  private

  public :: Particle, Species, SpeciesTable
  public :: particle_generator, generate_particle_arrays
  public :: load_species_from_file, species_index_of
  public :: read_particle_conditions
  public :: build_particles_from_species_domain

  !------------------------------
  ! Macro-particle storage (SoA)
  !------------------------------
  type :: Particle
    real(real64), allocatable :: position1(:), position2(:), position3(:)
    real(real64), allocatable :: velocity1(:), velocity2(:), velocity3(:)
    real(real64)  :: charge             = 0.0_real64      !! Coulomb
    real(real64)  :: mass               = 0.0_real64      !! kg
    real(real64)  :: weight             = 1.0_real64      !! macro weight
    real(real64)  :: temp_eV            = 0.0_real64      !! eV
    real(real64)  :: fractional_density = 0.0_real64      !! init weight
    character(len=16) :: name = ''                         !! "[e]" "[H-]" etc.
    integer(int8)     :: index = 0                         !! species code
  end type Particle

  !------------------------------
  ! Species metadata for chemistry
  !------------------------------
  type :: Species
    character(len=:), allocatable :: name   !! label like "[e]" "[H-]" "[H]"
    real(real64) :: mass          = 0.0_real64   !! kg
    real(real64) :: charge        = 0.0_real64   !! Coulomb
    real(real64) :: Ti_eV         = 0.0_real64   !! temperature
    real(real64) :: frac_density  = 0.0_real64   !! weight for init
  end type Species

  type :: SpeciesTable
    type(Species), allocatable :: species(:)   !! 1..number_species
    integer :: number_species = 0
    integer :: tag_neg  = 0     !! [H-] if present
    integer :: tag_neu  = 0     !! [H]  if present
    integer :: tag_beam = 0     !! [eb] if present
  end type SpeciesTable

  interface Particle
    module procedure particle_generator
  end interface

  contains

  !------------------------------------------------------------------
  ! Constructor: empty Particle with metadata set
  !------------------------------------------------------------------
  type(Particle) function particle_generator(charge, mass, weight, temp_eV, fractional_density, name, index) result(p)
    real(real64), intent(in) :: charge, mass, weight, temp_eV, fractional_density
    character(len=*), intent(in) :: name
    integer(int8),    intent(in) :: index
    allocate(p%position1(0), p%position2(0), p%position3(0))
    allocate(p%velocity1(0), p%velocity2(0), p%velocity3(0))
    p%charge             = charge
    p%mass               = mass
    p%weight             = weight
    p%temp_eV            = temp_eV
    p%fractional_density = fractional_density
    p%name  = adjustl(name(1:min(len_trim(name), len(p%name))))
    p%index = index
  end function particle_generator

  !------------------------------------------------------------------
  ! Allocate/resize SoA arrays for N macro-particles (zero-initialized).
  ! Uses int64 for N to support very large counts.
  !------------------------------------------------------------------
  
  subroutine build_particles_from_species_domain(table_species, dom, parts)
        use mod_geometry, only: Domain
        type(SpeciesTable), intent(in)              :: table_species
        type(Domain),       intent(in)              :: dom
        type(Particle),     allocatable, intent(out):: parts(:)

        integer(int32) :: nb_ppc
        real(real64)   :: Te_eV, ne_m3, ngas_m3
        integer        :: k
        integer(int64) :: N
        real(real64)   :: frac, ncell_prod_r8

        if (table_species%number_species <= 0) error stop 'build_particles_from_species_domain: empty SpeciesTable'
        allocate(parts(table_species%number_species))

        call read_particle_conditions(find_conditions_file(), Te_eV, ne_m3, ngas_m3, nb_ppc)

        ncell_prod_r8 = real(dom%ncell_x1,real64) * real(dom%ncell_x2,real64) * real(dom%ncell_x3,real64)

        do k = 1, table_species%number_species
            frac = max(0.0_real64, table_species%species(k)%frac_density)
            if (nb_ppc < 0) then
            N = 0_int64
            else
            N = nint( real(abs(nb_ppc),real64) * ncell_prod_r8 * frac, int64 )
            end if

            parts(k)%charge             = table_species%species(k)%charge
            parts(k)%mass               = table_species%species(k)%mass
            parts(k)%weight             = 1.0_real64
            parts(k)%temp_eV            = table_species%species(k)%Ti_eV
            parts(k)%fractional_density = table_species%species(k)%frac_density
            parts(k)%name               = table_species%species(k)%name
            parts(k)%index              = int(min(k,127), int8)

            call generate_particle_arrays(parts(k), N)
        end do
    end subroutine build_particles_from_species_domain
  
    subroutine generate_particle_arrays(p, N)
        type(Particle), intent(inout) :: p
        integer(int64), intent(in)    :: N
        integer :: Nd                  ! default INTEGER
        Nd = int(N)                    ! safe if you’ve already range-checked N

        if (N < 0_int64) error stop 'generate_particle_arrays: N must be >= 0'

        if (allocated(p%position1)) then
            if (size(p%position1, kind=int64) /= N) then
            deallocate(p%position1, p%position2, p%position3, &
                        p%velocity1, p%velocity2, p%velocity3)
            end if
        end if
        if (.not. allocated(p%position1)) then
            allocate(p%position1(Nd), p%position2(Nd), p%position3(Nd))
            allocate(p%velocity1(Nd), p%velocity2(Nd), p%velocity3(Nd))
        end if
        if (N > 0_int64) then
            p%position1 = 0.0_real64; p%position2 = 0.0_real64; p%position3 = 0.0_real64
            p%velocity1 = 0.0_real64; p%velocity2 = 0.0_real64; p%velocity3 = 0.0_real64
        end if
    end subroutine generate_particle_arrays

  !------------------------------------------------------------------
  ! Append one species into the table (grows allocatable array)
  !------------------------------------------------------------------
  subroutine append_species(table_species, label, mass_amu, charge_e, Ti_eV, ni0_in)
    type(SpeciesTable), intent(inout) :: table_species
    character(len=*),   intent(in)    :: label
    real(real64),       intent(in)    :: mass_amu, charge_e, Ti_eV, ni0_in
    type(Species), allocatable :: tmp(:)
    integer :: n, oldn, idx

    if (.not. allocated(table_species%species)) then
      allocate(table_species%species(1))
      oldn = 0
    else
      oldn = size(table_species%species)
      allocate(tmp(oldn)); tmp = table_species%species
      deallocate(table_species%species)
      allocate(table_species%species(oldn+1))
      table_species%species(1:oldn) = tmp
      deallocate(tmp)
    end if

    n   = oldn + 1
    idx = n
    table_species%number_species = n

    if (allocated(table_species%species(idx)%name)) deallocate(table_species%species(idx)%name)
    allocate(character(len=len_trim(label)) :: table_species%species(idx)%name)
    table_species%species(idx)%name = trim(label)

    table_species%species(idx)%charge       = charge_e * elementary_charge
    table_species%species(idx)%mass         = mass_amu  * amu
    table_species%species(idx)%Ti_eV        = Ti_eV
    table_species%species(idx)%frac_density = max(0.0_real64, ni0_in)

    select case (trim(table_species%species(idx)%name))
    case ('[H-]'); table_species%tag_neg  = idx
    case ('[H]');  table_species%tag_neu  = idx
    case ('[eb]'); table_species%tag_beam = idx
    end select
  end subroutine append_species

  !------------------------------------------------------------------
  ! Load species (IONS/NEUT sections) from chemistry file
  !------------------------------------------------------------------
  subroutine load_species_from_file(table_species)
    type(SpeciesTable), intent(inout) :: table_species
    character(len=256)  :: filename_particles
    integer             :: fileIndex, ios_in
    integer             :: ierr_mpi, rank_local
    character(len=256)  :: kind_electric_particle
    character(len=1024) :: line
    character(len=256)  :: tok(16)
    integer             :: ntok
    character(len=256)  :: label
    real(real64)        :: mass_in, charge_in, Ti_eV_in, ni0_in
    integer             :: outUnit, ios_out
    logical             :: out_is_open
    character(len=*), parameter :: outFile_particles = './Outputs/particles.out'
    logical :: in_neutrals_section, neutrals_header_written

    call MPI_Comm_rank(MPI_COMM_WORLD, rank_local, ierr_mpi)

    filename_particles = trim(find_particles_file())

    if (allocated(table_species%species)) deallocate(table_species%species)
    table_species%number_species = 0
    table_species%tag_neg = 0; table_species%tag_neu = 0; table_species%tag_beam = 0

    out_is_open = .false.; in_neutrals_section = .false.; neutrals_header_written = .false.
    if (rank_local == 0) then
      call execute_command_line('mkdir -p ./Outputs', exitstat=ios_out)
      open(newunit=outUnit, file=outFile_particles, status='replace', action='write', iostat=ios_out)
      if (ios_out == 0) out_is_open = .true.
    end if

    if (rank_local == 0 .and. out_is_open) then
      write(outUnit,'(A)') 'Electron:'
    end if
    call append_species(table_species, '[e]', mass_electron/amu, -1.0_real64, 5.0_real64, 1.0_real64)
    if (rank_local == 0 .and. out_is_open) then
      write(outUnit,'(A)') '[e]'
      write(outUnit,'(A)') ' '
      write(outUnit,'(A)') 'Ions:'
    end if

    open(newunit=fileIndex, file=filename_particles, status='old', action='read', iostat=ios_in)
    if (ios_in /= 0) then
      if (rank_local == 0) write(*,*) 'Could not open species file: ', trim(filename_particles)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr_mpi)
    end if

    do
      read(fileIndex,*,iostat=ios_in) kind_electric_particle
      if (ios_in /= 0) exit

      select case (kind_electric_particle(1:4))
      case ('NEUT','Neut','neut'); in_neutrals_section = .true.
      case ('IONS','Ions','ions'); in_neutrals_section = .false.
      case default; cycle
      end select

      do
        read(fileIndex,'(A)',iostat=ios_in) line
        if (ios_in /= 0) exit
        if (index(line,'----') > 0) exit
      end do
      if (ios_in /= 0) exit

      do
        read(fileIndex,'(A)',iostat=ios_in) line
        if (ios_in /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (index(line,'----') > 0) exit
        if (len_trim(line) >= 2 .and. line(1:2) == '--') cycle

        call split_ws(line, tok, ntok)
        if (ntok < 5) cycle

        label = tok(1)
        read(tok(2),*,iostat=ios_in) mass_in;   if (ios_in /= 0) cycle
        read(tok(3),*,iostat=ios_in) charge_in; if (ios_in /= 0) cycle
        read(tok(4),*,iostat=ios_in) Ti_eV_in;  if (ios_in /= 0) cycle
        read(tok(5),*,iostat=ios_in) ni0_in;    if (ios_in /= 0) cycle

        call append_species(table_species, trim(label), mass_in, charge_in, Ti_eV_in, ni0_in)

        if (rank_local == 0 .and. out_is_open) then
          if (in_neutrals_section .and. .not. neutrals_header_written) then
            write(outUnit,'(A)') ' '
            write(outUnit,'(A)') 'Neutrals:'
            neutrals_header_written = .true.
          end if
          write(outUnit,'(A)') trim(table_species%species(table_species%number_species)%name)
        end if
      end do
    end do

    close(fileIndex)

    if (rank_local == 0) then
      write(*,*) 'Loaded species: ', table_species%number_species
      if (table_species%tag_neg  > 0) write(*,*) 'tag_neg  -> ', trim(table_species%species(table_species%tag_neg)%name)
      if (table_species%tag_neu  > 0) write(*,*) 'tag_neu  -> ', trim(table_species%species(table_species%tag_neu)%name)
      if (table_species%tag_beam > 0) write(*,*) 'tag_beam -> ', trim(table_species%species(table_species%tag_beam)%name)
    end if
  end subroutine load_species_from_file

  !----------------------
  ! Small helpers
  !----------------------
  pure integer function species_index_of(table_species, label) result(idx)
    type(SpeciesTable), intent(in) :: table_species
    character(len=*),   intent(in) :: label
    integer :: k
    idx = 0
    do k=1,table_species%number_species
      if (trim(table_species%species(k)%name) == trim(label)) then
        idx = k; return
      end if
    end do
  end function species_index_of

  !------------------------------------------------------------------
  ! Read basic particle-related conditions from a small text file.
  !------------------------------------------------------------------
  subroutine read_particle_conditions(filename, electronTemperature, electronDensity, gasDensity, particlePerCell)
    character(*),   intent(in)  :: filename
    real(real64),   intent(out) :: electronTemperature, electronDensity, gasDensity
    integer(int32), intent(out) :: particlePerCell
    integer :: fileIndex, ios
    character(len=512) :: dummy

    electronTemperature = 0.0_real64
    electronDensity     = 0.0_real64
    gasDensity          = 0.0_real64
    particlePerCell     = 0_int32

    open(newunit=fileIndex, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'read_particle_conditions: cannot open file: '//trim(filename)

    read(fileIndex,'(A)', iostat=ios) dummy
    if (ios /= 0) error stop 'read_particle_conditions: unexpected EOF before line 2.'

    read(fileIndex,*,iostat=ios) electronTemperature
    if (ios /= 0) error stop 'read_particle_conditions: failed ElectronTemperature (line 2).'
    read(fileIndex,*,iostat=ios) electronDensity
    if (ios /= 0) error stop 'read_particle_conditions: failed ElectronDensity (line 3).'
    read(fileIndex,*,iostat=ios) gasDensity
    if (ios /= 0) error stop 'read_particle_conditions: failed GasDensity (line 4).'
    read(fileIndex,*,iostat=ios) particlePerCell
    if (ios /= 0) error stop 'read_particle_conditions: failed np_cell (line 5).'

    close(fileIndex)
  end subroutine read_particle_conditions

end module mod_particles
