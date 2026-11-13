module grid_monochromatic

  ! The purpose of this module is to take care of the emission of photons for
  ! the monochromatic radiative transfer code. The
  ! setup_monochromatic_grid_pdfs subroutine takes care of pre-computing the
  ! probability distribution function for emission at the frequency index inu
  ! (using the frequencies from the radiative transfer settings), and by
  ! default the emission is weighted by the total luminosity from each cell at
  ! that frequency. In future it will be easy to add an option to weigh things
  ! differently, for example by total energy in the cell, or giving equal
  ! weight to each cell. The emit_from_monochromatic_grid_pdf function is then
  ! used to emit a photon from the pre-computed probability distribution
  ! function. The allocate_ and deallocate_ subroutines should be called at
  ! the start and end of the monochormatic radiative transfer respectively.

  use core_lib
  use type_photon, only : photon
  use grid_geometry, only : geo, grid_sample_pdf_map, random_position_cell
  use type_dust, only : dust_sample_emit_probability, dust_emit
  use dust_main, only : d, n_dust, prepare_photon, update_optconsts
  use grid_physics
  use settings, only : frequencies

  implicit none
  save

  private
  public :: allocate_monochromatic_grid_pdfs
  public :: deallocate_monochromatic_grid_pdfs
  public :: setup_monochromatic_grid_pdfs
  public :: emit_from_monochromatic_grid_pdf
  public :: setup_dust_photon_pdfs
  public :: emit_dust_photon

  ! Variables for the monochromatic grid emission
  integer :: inu_current
  real(dp), allocatable :: mean_prob(:)
  type(pdf_discrete_dp), allocatable :: emiss_pdf(:)

  ! Polychromatic dust-photon pass (HYPERION_DUST_PHOTONS): cells sampled in
  ! proportion to their emitted luminosity, dust types in proportion to theirs
  type(pdf_discrete_dp), allocatable :: lum_pdf(:)
  real(dp), allocatable :: lum_type(:)
  ! optional biased component (HYPERION_DUST_PHOTONS_BIAS): cells in
  ! proportion to their luminosity shortward of a cut wavelength
  type(pdf_discrete_dp), allocatable :: lumb_pdf(:)
  real(dp), allocatable :: lumb_type(:), lum_cell(:,:), lumb_cell(:,:)
  real(dp) :: bias_frac = 0._dp

contains

  subroutine allocate_monochromatic_grid_pdfs()
    implicit none
    allocate(emiss_pdf(n_dust), mean_prob(n_dust))
  end subroutine allocate_monochromatic_grid_pdfs

  subroutine deallocate_monochromatic_grid_pdfs()
    implicit none
    deallocate(emiss_pdf, mean_prob)
  end subroutine deallocate_monochromatic_grid_pdfs

  subroutine setup_monochromatic_grid_pdfs(inu, empty)

    ! Sets up the probability distribution functions to emit from for a specific
    ! frequency for all dust types. This subroutine takes nu, the frequency to
    ! compute the PDF at, and modifies local variables in the module

    implicit none

    integer, intent(in) :: inu

    real(dp) :: nu
    integer :: dust_id
    integer :: emiss_var_id
    real(dp) :: emiss_var_frac
    real(dp), allocatable :: energy(:), prob(:)
    integer :: icell
    logical,intent(out) :: empty

    nu = frequencies(inu)

    ! Allocate temporary arrays
    allocate(energy(geo%n_cells), prob(geo%n_cells))

    ! Loop over dust types
    do dust_id=1, n_dust

       ! Find the total energy emitted inside each cell
       if (energy_abs_tot(dust_id) > 0._dp) then
          energy = specific_energy(:, dust_id) &
               &   * density(:, dust_id) &
               &   * geo%volume(:) &
               &   * dble(geo%n_cells) &
               &   / energy_abs_tot(dust_id)
       else
          energy = 0._dp
       end if

       ! Ensure that energy is zero in masked cells. Note that we can just
       ! work with all the cells here because the masked cells will get
       ! dropped out of the PDF anyway.
       if(geo%masked) then
          where(.not.geo%mask)
             energy = 0._dp
          end where
       end if

       ! Find in each cell the probability of emission at the desired frequency
       ! from the normalized emissivity PDF.
       do icell=1, geo%n_cells
          emiss_var_id = jnu_var_id(icell, dust_id)
          emiss_var_frac = jnu_var_frac(icell, dust_id)
          call dust_sample_emit_probability(d(dust_id), &
               &                            emiss_var_id, emiss_var_frac, &
               &                            nu, prob(icell))
       end do

       ! Set the PDF to the probability of emission at the frequency nu
       mean_prob(dust_id) = mean(prob * energy)
       if(mean_prob(dust_id) > 0._dp) call set_pdf(emiss_pdf(dust_id), prob * energy)

    end do

    ! Deallocate temporary arrays
    deallocate(energy, prob)

    inu_current = inu

    empty = sum(mean_prob) == 0._dp

  end subroutine setup_monochromatic_grid_pdfs

  type(photon) function emit_from_monochromatic_grid_pdf(inu) result(p)

    ! Given a frequency index (inu), sample a position in the grid using the
    ! locally pre-computed PDFs and set up a photon. This function is meant to
    ! be used in conjunction with setup_monochromatic_grid_pdfs, which, given
    ! a frequency, pre-computes the probabilty distribution functon for
    ! emission from the grid.

    implicit none

    integer, intent(in) :: inu

    integer :: dust_id
    real(dp) :: xi

    if(inu /= inu_current) call error('emit_from_monochromatic_grid_pdf', 'incorrect inu')

    p%nu = frequencies(inu)
    p%inu = inu

    call prepare_photon(p)
    call update_optconsts(p)

    call random(xi)
    dust_id = max(ceiling(xi*real(n_dust, dp)), 1)

    if(mean_prob(dust_id) == 0._dp) then
        p%energy = 0._dp
        return
    end if

    call grid_sample_pdf_map(emiss_pdf(dust_id), p%icell)
    p%in_cell = .true.

    ! Find random position inside cell
    call random_position_cell(p%icell, p%r) ! can probably make this a function

    ! Sample an isotropic direction
    call random_sphere_angle3d(p%a)
    call angle3d_to_vector3d(p%a, p%v)

    ! Set stokes parameters to unpolarized light
    p%s = stokes_dp(1._dp, 0._dp, 0._dp, 0._dp)

    ! Set the energy to 1., and it will be scaled in the main routine
    p%energy = mean_prob(dust_id)

    p%scattered = .false.
    p%reprocessed = .true.
    p%last_isotropic = .true.
    p%dust_id = dust_id
    p%last = 'de'

  end function emit_from_monochromatic_grid_pdf

  real(dp) function frac_above(pdf, x0) result(f)
    ! Fraction of a continuous PDF (x ascending) above x0, from its CDF
    implicit none
    type(pdf_dp), intent(in) :: pdf
    real(dp), intent(in) :: x0
    integer :: i
    real(dp) :: c
    if(x0 <= pdf%x(1)) then
       f = 1._dp
    else if(x0 >= pdf%x(pdf%n)) then
       f = 0._dp
    else
       i = 1
       do while(pdf%x(i+1) <= x0)
          i = i + 1
       end do
       c = pdf%cdf(i) + (pdf%cdf(i+1) - pdf%cdf(i)) * (x0 - pdf%x(i)) / (pdf%x(i+1) - pdf%x(i))
       f = max(0._dp, 1._dp - c / pdf%cdf(pdf%n))
    end if
  end function frac_above

  real(dp) function setup_dust_photon_pdfs(bias, bias_wav) result(lum_total)

    ! Luminosity-weighted emission PDFs for the polychromatic dust-photon
    ! pass of the final iteration: cell luminosity specific_energy * mass for
    ! each dust type (zero in masked cells and in cells at the minimum
    ! specific energy).  Returns the total luminosity
    ! (sum over dust types), which the caller divides among the photons.
    !
    ! bias > 0 (HYPERION_DUST_PHOTONS_BIAS): defensive-mixture importance
    ! sampling.  A fraction bias of the photons is drawn from cells in
    ! proportion to the luminosity they emit shortward of bias_wav microns
    ! (hot dust, the source of scattered near/mid-IR dust light), the rest in
    ! proportion to total luminosity; each photon carries the weight
    ! p_lum / p_mix (emit_dust_photon), so the estimate stays unbiased.

    implicit none

    real(dp), intent(in) :: bias, bias_wav
    real(dp), allocatable :: lum(:), fj(:)
    real(dp) :: nu_cut
    integer :: dust_id, k, ic, id

    if(allocated(lum_pdf)) deallocate(lum_pdf, lum_type)
    if(allocated(lumb_pdf)) deallocate(lumb_pdf, lumb_type)
    if(allocated(lum_cell)) deallocate(lum_cell, lumb_cell)
    allocate(lum_pdf(n_dust), lum_type(n_dust), lum(geo%n_cells))
    allocate(lumb_pdf(n_dust), lumb_type(n_dust))
    allocate(lum_cell(geo%n_cells, n_dust), lumb_cell(geo%n_cells, n_dust))
    lumb_type = 0._dp
    lumb_cell = 0._dp
    bias_frac = bias
    nu_cut = c_cgs / (bias_wav * 1.e-4_dp)

    do dust_id=1, n_dust
       lum = specific_energy(:, dust_id) * density(:, dust_id) * geo%volume(:)
       if(geo%masked) then
          where(.not.geo%mask)
             lum = 0._dp
          end where
       end if
       ! Cells held at the minimum specific energy emit energy they did not
       ! absorb (e.g. a 10 K ambient medium out to 1e6 AU); the absorb/re-emit
       ! chains of a stock final iteration never carry it, so leave it out
       ! (the raytracing pass still includes its direct emission)
       where(lum < 0._dp .or. specific_energy(:, dust_id) <= minimum_specific_energy(dust_id) * (1._dp + 1.e-10_dp))
          lum = 0._dp
       end where
       lum_type(dust_id) = sum(lum)
       if(lum_type(dust_id) > 0._dp) call set_pdf(lum_pdf(dust_id), lum)
       lum_cell(:, dust_id) = lum
       if(bias > 0._dp .and. lum_type(dust_id) > 0._dp) then
          allocate(fj(d(dust_id)%n_jnu))
          do k=1, d(dust_id)%n_jnu
             fj(k) = frac_above(d(dust_id)%j_nu(k), nu_cut)
          end do
          do ic=1, geo%n_cells
             if(lum(ic) > 0._dp) then
                id = jnu_var_id(ic, dust_id)
                lumb_cell(ic, dust_id) = lum(ic) * max(0._dp, fj(id) + jnu_var_frac(ic, dust_id) * (fj(id+1) - fj(id)))
             end if
          end do
          deallocate(fj)
          lumb_type(dust_id) = sum(lumb_cell(:, dust_id))
          if(lumb_type(dust_id) > 0._dp) call set_pdf(lumb_pdf(dust_id), lumb_cell(:, dust_id))
       end if
    end do

    deallocate(lum)

    lum_total = sum(lum_type)
    if(sum(lumb_type) <= 0._dp) bias_frac = 0._dp

  end function setup_dust_photon_pdfs

  type(photon) function emit_dust_photon() result(p)

    ! Emit a polychromatic dust photon: dust type and cell drawn in proportion
    ! to luminosity, frequency from the cell emissivity, isotropic direction.
    ! The energy is set to 1 and scaled by the caller.

    implicit none

    integer :: dust_id, id
    real(dp) :: xi, energy_scaling, acc, p_lum, p_mix
    logical :: biased

    ! mixture component (only with HYPERION_DUST_PHOTONS_BIAS, so that the
    ! unbiased pass draws the same random numbers as before)
    biased = .false.
    if(bias_frac > 0._dp) then
       call random(xi)
       biased = xi < bias_frac
    end if

    ! dust type in proportion to its luminosity (types that emit nothing have
    ! no PDF and are never chosen)
    call random(xi)
    dust_id = 0
    acc = 0._dp
    if(biased) then
       xi = xi * sum(lumb_type)
       do id=1, n_dust
          if(lumb_type(id) > 0._dp) then
             dust_id = id
             acc = acc + lumb_type(id)
             if(xi < acc) exit
          end if
       end do
    else
       xi = xi * sum(lum_type)
       do id=1, n_dust
          if(lum_type(id) > 0._dp) then
             dust_id = id
             acc = acc + lum_type(id)
             if(xi < acc) exit
          end if
       end do
    end if
    if(dust_id == 0) call error('emit_dust_photon', 'no dust luminosity')

    if(biased) then
       call grid_sample_pdf_map(lumb_pdf(dust_id), p%icell)
    else
       call grid_sample_pdf_map(lum_pdf(dust_id), p%icell)
    end if
    p%in_cell = .true.
    call random_position_cell(p%icell, p%r)

    p%dust_id = dust_id
    p%emiss_type = 3
    p%emiss_var_id = jnu_var_id(p%icell%ic, dust_id)
    p%emiss_var_frac = jnu_var_frac(p%icell%ic, dust_id)

    call dust_emit(d(dust_id), p%emiss_var_id, p%emiss_var_frac, p%nu, p%a, p%s, energy_scaling)
    call angle3d_to_vector3d(p%a, p%v)

    call prepare_photon(p)
    call update_optconsts(p)

    ! relative weight p_lum / p_mix (1 without bias); the caller multiplies by
    ! the total luminosity over the number of photons
    if(bias_frac > 0._dp) then
       p_lum = lum_cell(p%icell%ic, dust_id) / sum(lum_type)
       p_mix = (1._dp - bias_frac) * p_lum + bias_frac * lumb_cell(p%icell%ic, dust_id) / sum(lumb_type)
       p%energy = p_lum / p_mix
    else
       p%energy = 1._dp
    end if

    p%scattered = .false.
    p%reprocessed = .true.
    p%last_isotropic = .true.
    p%last = 'de'

  end function emit_dust_photon

end module grid_monochromatic
