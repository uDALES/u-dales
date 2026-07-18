program test_calc_psi
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ny = 1023, nz = 1023
    integer, parameter :: max_ly = 128, max_lz = 128
    integer, parameter :: number_of_repeats = 5

    integer :: iy, iz, j, k, nry, nrz, unit_no, repeat
    integer(int64) :: clock_start, clock_finish, clock_rate
    integer :: nly2(0:nz), nlz2(0:nz)
    real(dp) :: norm_y, norm_z
    real(dp) :: max_error, rms_error, relative_l2, tolerance
    real(dp) :: direct_time, fft_time, direct_average, fft_average
    real(dp) :: speedup, direct_checksum, fft_checksum
    real(dp), allocatable :: by(:,:), bz(:,:), randn(:,:)
    real(dp), allocatable :: psi_direct(:,:), psi_fft(:,:), error(:,:)

    nry = ny + 2*max(max_ly, max_lz)
    nrz = nz + 2*max(max_ly, max_lz)

    allocate(by(0:max_ly,0:nz), bz(0:max_lz,0:nz))
    allocate(randn(0:nry,0:nrz))
    allocate(psi_direct(0:ny,0:nz), psi_fft(0:ny,0:nz))
    allocate(error(0:ny,0:nz))

    ! Different filter lengths at different z levels exercise the same
    ! z-dependent bounds used by the production calc_psi routine.
    do iz = 0, nz
        nly2(iz) = max_ly ! 2 + mod(iz, max_ly - 1)
        nlz2(iz) = max_lz !1 + mod(iz, max_lz)
    end do

    ! Construct normalized, exponentially decaying filter coefficients.
    by = 0.0_dp
    bz = 0.0_dp
    do iz = 0, nz
        do j = 0, nly2(iz)
            by(j,iz) = exp(-real(j,dp)/(1.3_dp + 0.07_dp*real(iz,dp)))
        end do
        norm_y = sqrt(sum(by(0:nly2(iz),iz)**2))
        by(0:nly2(iz),iz) = by(0:nly2(iz),iz)/norm_y

        do k = 0, nlz2(iz)
            bz(k,iz) = exp(-real(k,dp)/(1.1_dp + 0.05_dp*real(iz,dp)))
        end do
        norm_z = sqrt(sum(bz(0:nlz2(iz),iz)**2))
        bz(0:nlz2(iz),iz) = bz(0:nlz2(iz),iz)/norm_z
    end do

    ! A deterministic field makes every run exactly reproducible.  Both
    ! methods below receive this same array.
    do iz = 0, nrz
        do iy = 0, nry
            randn(iy,iz) = sin(0.173_dp*real(iy + 1,dp) + &
                                  0.117_dp*real(iz + 1,dp)) + &
                            cos(0.071_dp*real(2*iy - iz,dp))
        end do
    end do

    ! Repeat each method so this small test takes long enough to time reliably.
    ! The checksums ensure that every repeated result remains observable.
    direct_checksum = 0.0_dp
    call system_clock(clock_start, clock_rate)
    do repeat = 1, number_of_repeats
        call calc_psi_direct(nly2, nlz2, by, bz, randn, psi_direct)
        direct_checksum = direct_checksum + &
            psi_direct(mod(7*repeat,ny+1),mod(3*repeat,nz+1))
    end do
    call system_clock(clock_finish)
    direct_time = real(clock_finish-clock_start,dp)/real(clock_rate,dp)

    fft_checksum = 0.0_dp
    call system_clock(clock_start, clock_rate)
    do repeat = 1, number_of_repeats
        call calc_psi_fft(nly2, nlz2, by, bz, randn, psi_fft)
        fft_checksum = fft_checksum + &
            psi_fft(mod(7*repeat,ny+1),mod(3*repeat,nz+1))
    end do
    call system_clock(clock_finish)
    fft_time = real(clock_finish-clock_start,dp)/real(clock_rate,dp)

    direct_average = direct_time/real(number_of_repeats,dp)
    fft_average = fft_time/real(number_of_repeats,dp)
    speedup = direct_time/max(fft_time,tiny(1.0_dp))

    error = abs(psi_fft - psi_direct)
    max_error = maxval(error)
    rms_error = sqrt(sum(error**2)/real(size(error),dp))
    relative_l2 = sqrt(sum((psi_fft - psi_direct)**2)) / &
                  max(sqrt(sum(psi_direct**2)), tiny(1.0_dp))
    tolerance = 1.0e-11_dp

    open(newunit=unit_no, file='calc_psi_comparison.csv', &
         status='replace', action='write')
    write(unit_no,'(A)') 'iy,iz,psi_direct,psi_fft,absolute_error'
    do iz = 0, nz
        do iy = 0, ny
            write(unit_no,'(I0,A,I0,3(A,ES24.16E3))') &
                iy, ',', iz, ',', psi_direct(iy,iz), ',', &
                psi_fft(iy,iz), ',', error(iy,iz)
        end do
    end do
    close(unit_no)

    write(*,'(A)')              'calc_psi direct-versus-FFT test'
    write(*,'(A,ES12.4E3)')     'maximum absolute error = ', max_error
    write(*,'(A,ES12.4E3)')     'RMS absolute error     = ', rms_error
    write(*,'(A,ES12.4E3)')     'relative L2 error      = ', relative_l2
    write(*,'(A,ES12.4E3)')     'acceptance tolerance   = ', tolerance
    write(*,'(A)') ''
    write(*,'(A,I0)')           'number of timed calls  = ', number_of_repeats
    write(*,'(A,F12.6,A)')      'direct total time      = ', direct_time, ' s'
    write(*,'(A,F12.6,A)')      'FFT total time         = ', fft_time, ' s'
    write(*,'(A,F12.6,A)')      'direct time per call   = ', direct_average, ' s'
    write(*,'(A,F12.6,A)')      'FFT time per call      = ', fft_average, ' s'
    write(*,'(A,F12.3,A)')      'FFT speed-up           = ', speedup, ' x'
    write(*,'(A,ES12.4E3)')     'direct checksum        = ', direct_checksum
    write(*,'(A,ES12.4E3)')     'FFT checksum           = ', fft_checksum

    if (max_error <= tolerance) then
        write(*,'(A)') 'PASS: calc_psi_FFT reproduces the direct sum.'
        write(*,'(A)') 'Details: calc_psi_comparison.csv'
    else
        write(*,'(A)') 'FAIL: the two calculations differ too much.'
        write(*,'(A)') 'Inspect calc_psi_comparison.csv.'
        stop 1
    end if

contains

    subroutine calc_psi_direct(nly2, nlz2, by, bz, randn, psi)
        implicit none

        integer, intent(in) :: nly2(0:), nlz2(0:)
        real(dp), intent(in) :: by(0:,0:), bz(0:,0:), randn(0:,0:)
        real(dp), intent(out) :: psi(0:,0:)
        integer :: iy, iz, j, k

        psi = 0.0_dp
        do iz = 0, ubound(psi,2)
            do iy = 0, ubound(psi,1)
                do k = 0, nlz2(iz)
                    do j = 0, nly2(iz)
                        psi(iy,iz) = psi(iy,iz) + &
                            by(j,iz)*bz(k,iz)*randn(iy+j,iz+k)
                    end do
                end do
            end do
        end do
    end subroutine calc_psi_direct


    ! subroutine calc_psi_fft(nly2, nlz2, by, bz, randn, psi)
        ! implicit none

        ! include 'fftw3.f'

        ! integer, parameter :: plan_kind = selected_int_kind(18)
        ! integer, intent(in) :: nly2(0:), nlz2(0:)
        ! real(dp), intent(in) :: by(0:,0:), bz(0:,0:), randn(0:,0:)
        ! real(dp), intent(out) :: psi(0:,0:)

        ! integer :: iy, iz, j, k, nfft, data_length, required_length
        ! integer(plan_kind) :: plan_data_forward, plan_kernel_forward
        ! integer(plan_kind) :: plan_data_backward
        ! complex(dp), allocatable :: data_hat(:), kernel_hat(:)

        ! data_length = size(randn,1)
        ! required_length = data_length + maxval(nly2)

        ! nfft = 1
        ! do while (nfft < required_length)
            ! nfft = 2*nfft
        ! end do

        ! allocate(data_hat(0:nfft-1), kernel_hat(0:nfft-1))

        ! call dfftw_plan_dft_1d(plan_data_forward, nfft, data_hat, &
                               ! data_hat, FFTW_FORWARD, FFTW_ESTIMATE)
        ! call dfftw_plan_dft_1d(plan_kernel_forward, nfft, kernel_hat, &
                               ! kernel_hat, FFTW_FORWARD, FFTW_ESTIMATE)
        ! call dfftw_plan_dft_1d(plan_data_backward, nfft, data_hat, &
                               ! data_hat, FFTW_BACKWARD, FFTW_ESTIMATE)

        ! psi = 0.0_dp

        ! do iz = 0, ubound(psi,2)
            ! ! Fourier transform of the y filter for this output z level.
            ! kernel_hat = cmplx(0.0_dp, 0.0_dp, kind=dp)
            ! do j = 0, nly2(iz)
                ! kernel_hat(j) = cmplx(by(j,iz), 0.0_dp, kind=dp)
            ! end do
            ! call dfftw_execute_dft(plan_kernel_forward, kernel_hat, kernel_hat)

            ! do k = 0, nlz2(iz)
                ! ! Transform one y row of the common random field.
                ! data_hat = cmplx(0.0_dp, 0.0_dp, kind=dp)
                ! do iy = 0, data_length - 1
                    ! data_hat(iy) = cmplx(randn(iy,iz+k), 0.0_dp, kind=dp)
                ! end do
                ! call dfftw_execute_dft(plan_data_forward, data_hat, data_hat)

                ! ! Multiplication by conjg(kernel_hat) implements
                ! ! sum_j by(j,iz)*randn(iy+j,iz+k), i.e. correlation.
                ! data_hat = data_hat*conjg(kernel_hat)
                ! call dfftw_execute_dft(plan_data_backward, data_hat, data_hat)

                ! do iy = 0, ubound(psi,1)
                    ! psi(iy,iz) = psi(iy,iz) + bz(k,iz)* &
                        ! real(data_hat(iy),dp)/real(nfft,dp)
                ! end do
            ! end do
        ! end do

        ! call dfftw_destroy_plan(plan_data_forward)
        ! call dfftw_destroy_plan(plan_kernel_forward)
        ! call dfftw_destroy_plan(plan_data_backward)
        ! deallocate(data_hat, kernel_hat)
    ! end subroutine calc_psi_fft
	
	subroutine calc_psi_fft(nly2, nlz2, by, bz, randn, psi)
		implicit none

		include 'fftw3.f'

		integer, parameter :: plan_kind = selected_int_kind(18)

		integer, intent(in) :: nly2(0:), nlz2(0:)
		real(dp), intent(in) :: by(0:,0:), bz(0:,0:)
		real(dp), intent(in) :: randn(0:,0:)
		real(dp), intent(out) :: psi(0:,0:)

		integer :: iy, iz, j, k, nfft
		integer, save :: cached_nfft = 0
		logical, save :: initialized = .false.

		integer(plan_kind), save :: plan_forward
		integer(plan_kind), save :: plan_backward

		real(dp), allocatable, save :: work_real(:)
		complex(dp), allocatable, save :: work_hat(:)
		complex(dp), allocatable, save :: kernel_hat(:)

		! The supplied y halo already contains every randn(iy+j,...)
		! required for the requested output range.
		nfft = size(randn,1)

		if (nfft < size(psi,1) + maxval(nly2)) then
			error stop 'Insufficient y halo in calc_psi_fft'
		end if

		! Cache allocation and FFTW plans between calls.
		if (.not. initialized .or. cached_nfft /= nfft) then
			if (initialized) then
				call dfftw_destroy_plan(plan_forward)
				call dfftw_destroy_plan(plan_backward)
				deallocate(work_real, work_hat, kernel_hat)
			end if

			allocate(work_real(0:nfft-1))
			allocate(work_hat(0:nfft/2))
			allocate(kernel_hat(0:nfft/2))

			call dfftw_plan_dft_r2c_1d(plan_forward, nfft,         &
									   work_real, work_hat,        &
									   FFTW_MEASURE)

			call dfftw_plan_dft_c2r_1d(plan_backward, nfft,        &
									   work_hat, work_real,        &
									   FFTW_MEASURE)

			cached_nfft = nfft
			initialized = .true.
		end if

		psi = 0.0_dp

		do iz = 0, ubound(psi,2)

			! Transform the z-dependent y filter.
			work_real = 0.0_dp
			do j = 0, nly2(iz)
				work_real(j) = by(j,iz)
			end do

			call dfftw_execute_dft_r2c(plan_forward,               &
									   work_real, work_hat)
			kernel_hat = work_hat

			! Use linearity to perform the complete z filtering before
			! the y-direction FFT:
			!
			! sum_k bz(k)*FFT(randn(:,iz+k))
			!   = FFT(sum_k bz(k)*randn(:,iz+k))
			work_real = 0.0_dp
			do k = 0, nlz2(iz)
				work_real = work_real + bz(k,iz) *                 &
							randn(0:nfft-1,iz+k)
			end do

			call dfftw_execute_dft_r2c(plan_forward,               &
									   work_real, work_hat)

			! conjg(kernel_hat) implements correlation in y:
			! sum_j by(j,iz)*randn(iy+j,...)
			work_hat = work_hat * conjg(kernel_hat)

			call dfftw_execute_dft_c2r(plan_backward,              &
									   work_hat, work_real)

			do iy = 0, ubound(psi,1)
				psi(iy,iz) = work_real(iy) / real(nfft,dp)
			end do
		end do

	end subroutine calc_psi_fft

end program test_calc_psi
