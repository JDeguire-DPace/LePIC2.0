module mod_particles
    !! =========================================================================
    !!  Module: mod_particles
    !!  - Particle (SoA container for macro-particles)
    !!  - Species / SpeciesTable metadata used by collisions
    !!  - particle_generator
    !!  - generate_particle_arrays
    !!  - load_species_from_file (parses IONS/NEUT sections)
    !!  - species_index_of
    !!  - read_particle_conditions
    !!  - build_particles_from_species  <--- NEW (allocates arrays per species)
    !! =========================================================================
    use iso_fortran_env, only: real64, int32, int8, int64
    use mpi
    use mod_InputFiles
    use mod_functionsText
    use mod_constants                ! expects: elementary_charge, amu, etc.
    implicit none
    private

    !------------------------------
    ! Public API
    !------------------------------
    public :: Particle, Species, SpeciesTable
    public :: particle_generator, generate_particle_arrays
    public :: load_species_from_file, species_index_of
    public :: read_particle_conditions
    public :: build_particles_from_species

    !------------------------------
    ! Macro-particle storage (SoA)
    !------------------------------
    type :: Particle
        real(real64), allocatable :: position1(:), position2(:), position3(:)
        real(real64), allocatable :: velocity1(:), velocity2(:), velocity3(:)
        real(real64)  :: charge = 0.0_real64      !! Coulomb
        real(real64)  :: mass   = 0.0_real64      !! kg
        real(real64)  :: weight = 1.0_real64      !! macro weight (particles per superparticle)
        real(real64)  :: temp_eV = 0.0_real64     !! eV
        real(real64)  :: fractional_density = 0.0_real64 !! species weight for init
        character(len=16) :: name = ''            !! e.g., "[e]", "[H-]"
        integer(int8)     :: index = 0            !! species code
    end type Particle

    !------------------------------
    ! Species metadata for chemistry
    !------------------------------
    type :: Species
        character(len=:), allocatable :: name    !! label like "[e]" "[H-]" "[H]"
        real(real64)     :: mass = 0.0    !! kg
        real(real64)     :: charge = 0.0  !! Coulomb
        real(real64)     :: Ti_eV = 0.0   !! temperature
        real(real64)     :: frac_density = 0.0  !! weight for init (see notes)
    end type Species

    type :: SpeciesTable
        type(Species), allocatable :: species(:)     !! 1..number_species
        integer :: number_species = 0
        integer :: tag_neg  = 0                !! [H-] if present
        integer :: tag_neu  = 0                !! [H]  if present
        integer :: tag_beam = 0                !! [eb] if present
    end type SpeciesTable

    interface Particle
        module procedure particle_generator
    end interface

contains

    !--------------------------------------------------------------------------
    ! Constructor-style function: make an empty Particle with metadata set.
    !--------------------------------------------------------------------------
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

    !--------------------------------------------------------------------------
    ! Allocate/resize SoA arrays for N macro-particles (zero-initialized).
    !--------------------------------------------------------------------------
    subroutine generate_particle_arrays(p, N)
        type(Particle), intent(inout) :: p
        integer,        intent(in)    :: N
        if (N < 0) error stop 'generate_particle_arrays: N must be >= 0'
        if (allocated(p%position1)) then
            if (size(p%position1) /= N) then
                deallocate(p%position1, p%position2, p%position3, &
                           p%velocity1, p%velocity2, p%velocity3)
            end if
        end if
        if (.not. allocated(p%position1)) then
            allocate(p%position1(N), p%position2(N), p%position3(N))
            allocate(p%velocity1(N), p%velocity2(N), p%velocity3(N))
        end if
        if (N > 0) then
            p%position1 = 0.0_real64; p%position2 = 0.0_real64; p%position3 = 0.0_real64
            p%velocity1 = 0.0_real64; p%velocity2 = 0.0_real64; p%velocity3 = 0.0_real64
        end if
    end subroutine generate_particle_arrays

    !--------------------------------------------------------------------------
    ! Load species (IONS/NEUT sections) from chemistry file into SpeciesTable.
    !--------------------------------------------------------------------------
    subroutine load_species_from_file(table_species)
        implicit none

        type(SpeciesTable), intent(inout) :: table_species

        ! --- local I/O and MPI ---
        character(len=256) :: filename_particles
        integer            :: fileIndex, ios_in
        integer            :: ierr_mpi, rank_local

        ! --- parsing helpers ---
        character(len=256) :: kind_electric_particle     ! "IONS" or "NEUT"
        character(len=1024):: line                       ! read whole line (long to be safe)
        character(len=256) :: tok(16)                    ! tokens after split (up to 16)
        integer            :: ntok
        character(len=256) :: label                      ! token(1) -> name (e.g., [DeuteriumV9])
        real(real64)       :: mass_in, charge_in, Ti_eV_in, ni0_in

        ! --- output file for species names (rank 0 only) ---
        integer                      :: outUnit
        logical                      :: out_is_open
        integer                      :: ios_out
        character(len=*), parameter  :: outFile_particles = './Outputs/particles.out'

        ! --- flags to insert a separator before neutrals ---
        logical :: in_neutrals_section       = .false.
        logical :: neutrals_header_written   = .false.
        ! =============================

        ! Get local MPI rank
        call MPI_Comm_rank(MPI_COMM_WORLD, rank_local, ierr_mpi)

        filename_particles = trim(find_particles_file())

        ! Start fresh
        if (allocated(table_species%species)) deallocate(table_species%species)
        table_species%number_species = 0
        table_species%tag_neg  = 0
        table_species%tag_neu  = 0
        table_species%tag_beam = 0

        ! Output prep
        out_is_open = .false.
        if (rank_local == 0) then
            call execute_command_line('mkdir -p ./Outputs', exitstat=ios_out)
            open(newunit=outUnit, file=outFile_particles, status='replace', action='write', iostat=ios_out)
            if (ios_out /= 0) then
                write(*,*) 'Warning: could not open ', trim(outFile_particles), ' (iostat=', ios_out, ').'
            else
                out_is_open = .true.
            end if
        end if

        ! Electron always first (optional)
        if (rank_local == 0 .and. out_is_open) then
            write(outUnit,'(A)') 'Electron:'
        end if
        call append_species(table_species, '[e]', mass_electron/amu, -1.0_real64, 5.0_real64, 1.0_real64)
        if (rank_local == 0 .and. out_is_open) then
            write(outUnit,'(A)') '[e]'
            write(outUnit,'(A)') ' '
            write(outUnit,'(A)') 'Ions:'
        end if

        ! Open species file
        open(newunit=fileIndex, file=trim(filename_particles), status='old', action='read', iostat=ios_in)
        if (ios_in /= 0) then
            if (rank_local == 0) write(*,*) 'Could not open species file: ', trim(filename_particles)
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr_mpi)
        end if

        ! Read species blocks
        do
            read(fileIndex,*,iostat=ios_in) kind_electric_particle
            if (ios_in /= 0) exit

            select case (kind_electric_particle(1:4))
            case ('NEUT','Neut','neut')
                in_neutrals_section = .true.
            case ('IONS','Ions','ions')
                in_neutrals_section = .false.
            case default
                cycle
            end select

            ! Skip until '----'
            do
                read(fileIndex,'(A)',iostat=ios_in) line
                if (ios_in /= 0) exit
                if (index(line,'----') > 0) exit
            end do
            if (ios_in /= 0) exit

            ! Read entries in this section
            do
                read(fileIndex,'(A)',iostat=ios_in) line
                if (ios_in /= 0) exit

                ! Trim leading spaces, ignore blanks, comments, or section breaks
                line = adjustl(line)
                if (len_trim(line) == 0) cycle
                if (index(line,'----') > 0) exit
                if (len_trim(line) >= 2) then
                    if (line(1:2) == '--') cycle
                end if

                ! Split by whitespace; expect at least 5 tokens
                call split_ws(line, tok, ntok)
                if (ntok < 5) cycle

                label = tok(1)
                read(tok(2),*,iostat=ios_in) mass_in
                if (ios_in /= 0) cycle
                read(tok(3),*,iostat=ios_in) charge_in
                if (ios_in /= 0) cycle
                read(tok(4),*,iostat=ios_in) Ti_eV_in
                if (ios_in /= 0) cycle
                read(tok(5),*,iostat=ios_in) ni0_in
                if (ios_in /= 0) cycle

                call append_species(table_species, trim(label), mass_in, charge_in, Ti_eV_in, ni0_in)

                ! Write to particles.out as we append
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

        ! Display summary on rank 0 only
        if (rank_local == 0) then
            write(*,*) 'Loaded species: ', table_species%number_species
            if (table_species%tag_neg  > 0) write(*,*) 'tag_neg  -> ', trim(table_species%species(table_species%tag_neg)%name)
            if (table_species%tag_neu  > 0) write(*,*) 'tag_neu  -> ', trim(table_species%species(table_species%tag_neu)%name)
            if (table_species%tag_beam > 0) write(*,*) 'tag_beam -> ', trim(table_species%species(table_species%tag_beam)%name)
        end if




    end subroutine load_species_from_file






    !======================
    ! Helpers
    !======================
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

    subroutine append_species(table_species, label, mass_amu, charge_e, Ti_eV, ni0_in)
        type(SpeciesTable), intent(inout) :: table_species
        character(len=*),   intent(in)    :: label
        real(real64),       intent(in)    :: mass_amu, charge_e, Ti_eV, ni0_in

        type(Species), allocatable :: tmp(:)
        integer :: n, oldn, idx

        ! grow array (simple push-back)
        if (.not. allocated(table_species%species)) then
            allocate(table_species%species(1))
            oldn = 0
        else
            oldn = size(table_species%species)
            allocate(tmp(oldn))
            tmp = table_species%species
            deallocate(table_species%species)
            allocate(table_species%species(oldn+1))
            table_species%species(1:oldn) = tmp
            deallocate(tmp)
        end if

        n   = oldn + 1
        idx = n
        table_species%number_species = n

        ! --- set name with exact length ---
        if (allocated(table_species%species(idx)%name)) deallocate(table_species%species(idx)%name)
        allocate(character(len=len_trim(label)) :: table_species%species(idx)%name)
        table_species%species(idx)%name = trim(label)

        ! --- set properties ---
        table_species%species(idx)%charge       = charge_e * elementary_charge
        table_species%species(idx)%mass         = mass_amu  * amu
        table_species%species(idx)%Ti_eV        = Ti_eV
        table_species%species(idx)%frac_density = max(0.0_real64, ni0_in)

        ! --- tags by name (adjust as you wish) ---
        select case (trim(table_species%species(idx)%name))
        case ('[H-]'); table_species%tag_neg  = idx
        case ('[H]');  table_species%tag_neu  = idx
        case ('[eb]'); table_species%tag_beam = idx
        end select
    end subroutine append_species


    !--------------------------------------------------------------------------
    ! Read basic particle-related conditions from a small text file.
    !--------------------------------------------------------------------------
    subroutine read_particle_conditions(filename, electronTemperature, electronDensity, gasDensity, particlePerCell)
        character(*), intent(in) :: filename
        real(real64), intent(out) :: electronTemperature, electronDensity, gasDensity
        integer(int32), intent(out) :: particlePerCell

        integer :: fileIndex, ios
        character(len=:), allocatable :: dummy

        electronTemperature = 0.0_real64
        electronDensity     = 0.0_real64
        gasDensity          = 0.0_real64
        particlePerCell     = 0_int32

        open(newunit=fileIndex, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'read_particle_conditions: cannot open file: '//trim(filename)

        ! Line 1: header/comment — discard
        read(fileIndex,'(A)', iostat=ios) dummy
        if (ios /= 0) error stop 'read_particle_conditions: unexpected EOF before line 2.'

        ! Lines 2–5: values (comments after "!" are ignored by list-directed IO)
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

    !--------------------------------------------------------------------------
    ! NEW: Build per-species Particle containers and allocate their SoA arrays.
    !   - table_species         : SpeciesTable filled by load_species_from_file
    !   - np_cell     : desired macro-particles per cell (baseline)
    !   - ncell_local : number of grid cells owned by THIS MPI rank
    !   - parts(:)    : output array of Particle, one per species
    !
    ! N_spec(k) = round( np_cell * ncell_local * frac_density(k) )
    !   * If frac_density=1 -> baseline np_cell per cell
    !   * If you use true fractions (e.g. 0.2/0.8), keep frac_density in [0,1]
    !--------------------------------------------------------------------------
    subroutine build_particles_from_species(table_species, parts)
        type(SpeciesTable), intent(in)           :: table_species
        type(Particle),     allocatable, intent(out) :: parts(:)

        integer(int32) :: nb_ppc
        real(real64)   :: temperature_electron
        real(real64)   :: initial_density
        real(real64)   :: gas_density

        ! >>> FIXED: these are integers, not integer(real64)
        integer(int32) :: ncell_x1, ncell_x2, ncell_x3
        ! <<<

        integer :: k
        integer(int64) :: N64
        integer        :: N
        real(real64)   :: frac

        ! --- read grid size from the first line of geometry file ---
        character(len=512) :: line, head
        character(len=512) :: geomfile
        integer :: u, ios, bang

        geomfile = trim(find_geometry_file())
        open(newunit=u, file=geomfile, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'build_particles_from_species: could not open geometry file: ', trim(geomfile), ' iostat=', ios
            error stop
        end if

        read(u,'(A)',iostat=ios) line
        close(u)
        if (ios /= 0) then
            write(*,*) 'build_particles_from_species: failed to read first line of ', trim(geomfile), ' iostat=', ios
            error stop
        end if

        ! strip trailing comment starting with '!'
        bang = index(line, '!')
        if (bang > 0) then
            head = line(:bang-1)
        else
            head = line
        end if

        ! read the first three integers (space-separated) from the head
        ncell_x1 = 0_int64; ncell_x2 = 0_int64; ncell_x3 = 0_int64
        read(head,*,iostat=ios) ncell_x1, ncell_x2, ncell_x3
        if (ios /= 0) then
            write(*,*) 'build_particles_from_species: could not parse ncell_x1/x2/x3 from: "', trim(head), '"'
            error stop
        end if
        ! now ncell_x1, ncell_x2, ncell_x3 are set from the geometry file

        if (table_species%number_species <= 0) error stop 'build_particles_from_species: empty SpeciesTable'

        allocate(parts(table_species%number_species))

        call read_particle_conditions(find_conditions_file(), temperature_electron, initial_density, gas_density, nb_ppc)

        do k = 1, table_species%number_species
            frac = max(0.0_real64, table_species%species(k)%frac_density)

            ! Compute count in 64-bit then convert safely to default INTEGER
            N64 = int( real(nb_ppc,real64) * real(ncell_x1*ncell_x2*ncell_x3,real64) * frac * initial_density, int64 )

            if (N64 > huge(N)) then
                write(*,*) 'Warning: particle count overflow; clamping for species ', trim(table_species%species(k)%name)
                N = huge(N)
            else
                N = int(N64)
            end if

            ! Copy metadata to runtime Particle container
            parts(k)%charge             = table_species%species(k)%charge
            parts(k)%mass               = table_species%species(k)%mass
            parts(k)%weight             = 1.0_real64
            parts(k)%temp_eV            = table_species%species(k)%Ti_eV
            parts(k)%fractional_density = table_species%species(k)%frac_density
            parts(k)%name               = table_species%species(k)%name
            parts(k)%index              = int(min(k, 127), int8)  ! safe int8 label

            call generate_particle_arrays(parts(k), N64)
        end do
    end subroutine build_particles_from_species


end module mod_particles