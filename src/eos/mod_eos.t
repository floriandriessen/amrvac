!> This module defines the handling of the eos of amrvac

module mod_eos
    use mod_global_parameters
    use mod_eos_container
    implicit none
    private

    character(len=std_len), parameter :: neOnh_table_prefix = "LTEeos_grid_neOnh_"

    !> eos object
    type(eos_container), public, allocatable     :: eos

    public :: eos_init

    contains
        !> Read this module"s parameters from a file
        subroutine eos_read_params(files)
            character(len=*), intent(in) :: files(:)
            integer                      :: n

            !> Default values
            character(len=std_len)  :: eos_type
            logical :: ionE
            logical :: table_check
            double precision :: He_abundance
            double precision :: gamma
            double precision :: gamma_minus_1
            double precision :: inv_gamma
            double precision :: inv_gamma_minus_1

            eos_type = 'FI'
            ionE = .false.
            table_check = .false.
            He_abundance = 0.1d0
            gamma = 5.0d0/3.0d0
            gamma_minus_1 = gamma - 1.0d0
            inv_gamma = 1.0d0 / gamma
            inv_gamma_minus_1 = 1.0d0 / gamma_minus_1

            namelist /eos_list/ eos_type, ionE, table_check, He_abundance, gamma

            do n = 1, size(files)
                open(unitpar, file=trim(files(n)), status="old")
                read(unitpar, eos_list, end=111)
        111     close(unitpar)
            end do

            eos%eos_type = eos_type
            eos%ionE = ionE
            if ((eos%eos_type == 'FI') .and. eos%ionE) then
                eos%ionE = .false.
                if(mype==0) write(*,*) "WARNING: eos_type = 'FI' requires ionE = .false."
            endif
            eos%table_check = table_check
            eos%He_abundance = He_abundance
            eos%gamma = gamma
            eos%gamma_minus_1 = eos%gamma - 1.0d0
            eos%inv_gamma = 1.0d0 / eos%gamma
            eos%inv_gamma_minus_1 = 1.0d0 / eos%gamma_minus_1
        
        end subroutine eos_read_params

        subroutine eos_init()
            
            allocate(eos)
            call eos_read_params(par_files)
            if (eos%eos_type == 'LTE') then
                call load_lte_tables()
            endif

        end subroutine eos_init

        !> This hasn't been tested yet - it needs an LTE use case, so far only considered Orszag_Tang
        !> This is deliberate as I want the FI case to work before I extend it.
        subroutine load_lte_tables()
            character(len=std_len) :: subname
            double precision, allocatable :: neOnH_table(:,:)
            
            subname = "H" !> Always consider at least hydrogen
            
            if (eos%He_abundance > 0.0d0) then
                write(subname,"(A,A)"), trim(subname), "HE"
            endif

            if (eos%ionE) then !> Consider the ionisation energy?
                write(subname,"(A,A)"), trim(subname), "_ionE"
            else
                write(subname,"(A,A)"), trim(subname), "_NoionE"
            endif

                
            write(eos%neOnH%filename,"(A,A,A)"), trim(neOnh_table_prefix), trim(subname), '.bin'
            if (mype==0) then
                print*, "Reading neOnH EoS table from: ", trim(eos%neOnH%filename)
            endif
            call read_eos_from_file(trim(eos%neOnH%filename), eos%neOnH%table, eos%neOnH%dim1, eos%neOnH%dim2, eos%neOnH%var1_min, eos%neOnH%var1_max, eos%neOnH%var2_min, eos%neOnH%var2_max, return_log=.true.) !> All variables provided in logspace

        end subroutine load_lte_tables

        subroutine initialise_eos_spaces()
            use mod_input_output

            if (eos%eos_type == 'LTE') then
                allocate(eos%T_space(ixG^T),eos%y_space(ixG^T),eos%p_to_eint_factor(ixG^T))
            else
                !> idea being we can just include this in the pre-existing routines
                !> might need to be careful about name so as to be readable
                allocate(eos%p_to_eint_factor(ixG^T)) !> actually equal to 1/(gamma-1) + ionE/pressure
            endif

            if (.not. eos%ionE) then
                eos%p_to_eint_factor(ixG^T) = eos%inv_gamma_minus_1 !> since ionE == 0
            endif

        end subroutine initialise_eos_spaces

        subroutine read_eos_from_file(filename, quantity_in, dimy, dimx, var1_min, var1_max, var2_min, var2_max, return_log)
            character(len=*), intent(in) :: filename
            double precision, intent(out) :: var1_min, var1_max, var2_min, var2_max
            integer, intent(out) :: dimx, dimy
            double precision, allocatable, intent(inout) :: quantity_in(:,:)
            logical, intent(in), optional :: return_log

            logical :: return_log_op

            if (present(return_log)) then
                return_log_op = return_log
            else
                return_log_op = .false.
            endif

            open(unit=10, file=filename, access='stream', form='unformatted', action='read')
            
            !> First read the header information, size from python file, and allocate
            read(10) dimy,dimx
            allocate(quantity_in(dimy,dimx))
            
            !> Then read the header information, bounds from python file
            read(10) var1_min, var1_max, var2_min, var2_max !> nhmin, nhmax, tmin, tmax
            
            !> Then read the 2D cube from the python file
            read(10) quantity_in

            if (mype==0) then
                if (eos%table_check) then
                    print*, filename, " table check"
                    print*, var1_min, var1_max, var2_min, var2_max
                    print*, 'corners', ' interpolated'
                    print*, dlog10(quantity_in(1,1)), interp_clamped_bilinear_table_log10(var1_min, var2_min, dlog10(quantity_in), dimy, dimx, var1_min, var1_max, var2_min, var2_max)
                    print*, dlog10(quantity_in(1,dimx)), interp_clamped_bilinear_table_log10(var1_min, var2_max, dlog10(quantity_in), dimy, dimx, var1_min, var1_max, var2_min, var2_max)
                    print*, dlog10(quantity_in(dimy,1)), interp_clamped_bilinear_table_log10(var1_max, var2_min, dlog10(quantity_in), dimy, dimx, var1_min, var1_max, var2_min, var2_max)
                    print*, dlog10(quantity_in(dimy,dimx)), interp_clamped_bilinear_table_log10(var1_max, var2_max, dlog10(quantity_in), dimy, dimx, var1_min, var1_max, var2_min, var2_max)
                    print*, "These should be identical, if not then the table structure is inconsistent"
                    print*, "###########"
                endif
            endif

            close(10)

            if (return_log_op) then
                quantity_in = dlog10(quantity_in)
            endif

        end subroutine read_eos_from_file

        pure double precision function interp_clamped_bilinear_table_log10(vary, varx, table, nx, ny, vmin_y, vmax_y, vmin_x, vmax_x) result(z)
            double precision, intent(in) :: vary, varx
            integer, intent(in)  :: nx, ny
            double precision, intent(in) :: table(ny, nx)
            double precision, intent(in) :: vmin_y, vmax_y, vmin_x, vmax_x

            double precision :: fx, fy, tx, ty, rx, ry, xstep_inv, ystep_inv
            integer  :: ix, iy, ix1, iy1
            
            xstep_inv = dble(nx-1) / (vmax_x - vmin_x)
            ystep_inv = dble(ny-1) / (vmax_y - vmin_y)

            ! fractional positions
            fx = (varx - vmin_x) * xstep_inv
            fy = (vary - vmin_y) * ystep_inv

            ! clamp to [0, nx-1] and [0, ny-1]
            rx = max(0.0d0, min(fx, dble(nx-1)))
            ry = max(0.0d0, min(fy, dble(ny-1)))

            ix = int(rx)               ! 0 .. nx-1
            iy = int(ry)               ! 0 .. ny-1

            ix1 = min(ix+1, nx-1)
            iy1 = min(iy+1, ny-1)

            tx = rx - dble(ix)      ! 0..1
            ty = ry - dble(iy)      ! 0..1

            z = (1.0d0-ty) * ((1.0d0-tx)*table(iy+1, ix+1) + tx*table(iy+1, ix1+1)) &
            +        ty  * ((1.0d0-tx)*table(iy1+1,ix+1) + tx*table(iy1+1,ix1+1))
        end function interp_clamped_bilinear_table_log10

end module mod_eos
!> Needs a line after to pass the preprocesor