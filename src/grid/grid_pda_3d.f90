! This approximation is not valid for cells which have emission in them, but
! that's fine, because if they had emission in them we wouldn't be calling this
! here

module grid_pda

  use core_lib
  use type_grid_cell
  use grid_physics
  use grid_geometry
  use dust_main
  use grid_io
  use grid_pda_geometry

  implicit none
  save

  private
  public :: solve_pda

  real(dp), parameter :: tolerance_iter = 1.e-4 ! energy calculation convergence criterion
  real(dp), parameter :: tolerance_exact = 1.e-5 ! energy calculation convergence criterion
  real(dp), parameter :: threshold_pda = 0.005 ! maximum number of photons required to use PDA

  real(dp), allocatable :: e_mean(:)

  ! HYPERION_PDA_SOLVER (read once): 'banded' solves the PDA system with a
  ! band LU factorisation for any number of PDA cells; any other value keeps
  ! the stock solvers and only adds a wall-time line per PDA solve; unset
  ! gives the stock behaviour and output.
  logical :: pda_env_read = .false.
  logical :: pda_banded = .false.
  logical :: pda_timed = .false.
  logical :: pda_report = .false.
  ! 'compare': every sweep also runs the banded solver on the same state,
  ! prints its difference from the stock solver, and continues with the stock
  ! result (validation)
  logical :: pda_compare = .false.
  ! largest band matrix (elements) the banded solver allocates; above it the
  ! stock solver runs instead
  real(dp), parameter :: pda_band_max_elements = 4.e8_dp

contains

  real(dp) elemental function difference_ratio(a, b)
    implicit none
    real(dp), intent(in) :: a, b
    difference_ratio = max(a/b, b/a)
  end function difference_ratio

  subroutine update_specific_energy(ic)
    implicit none
    integer,intent(in) :: ic
    integer :: id
    real(dp) :: s_prev, s, smin, smax
    do id=1,n_dust

       s = specific_energy(ic, id)

       smin = d(id)%specific_energy(1)
       smax = d(id)%specific_energy(d(id)%n_e)

       if(e_mean(ic) < smin / kappa_planck(id, smin)) then
          call warn("update_specific_energy", "specific energy in PDA below minimum allowed by dust type - resetting")
          s = smin
       else if (e_mean(ic) > smax / kappa_planck(id, smax)) then
          call warn("update_specific_energy", "specific energy in PDA above maximum allowed by dust type - resetting")
          s = smax
       else
          do
             s_prev = s
             s = e_mean(ic) * kappa_planck(id, s)
             if(difference_ratio(s, s_prev) - 1._dp < 1.e-5_dp) exit
          end do
       end if
       specific_energy(ic, id) = s
    end do
  end subroutine update_specific_energy

  subroutine update_e_mean(ic)
    implicit none
    integer,intent(in) :: ic
    integer :: id
    e_mean(ic) = 0.
    if(sum(density(ic, :)) > 0._dp) then
       do id=1,n_dust
          e_mean(ic) = e_mean(ic) + density(ic,id) * specific_energy(ic,id) / kappa_planck(id, specific_energy(ic,id))
       end do
       e_mean(ic) = e_mean(ic) / sum(density(ic, :))
    end if
  end subroutine update_e_mean

  subroutine solve_pda()

    implicit none

    real(dp),allocatable :: specific_energy_prev(:,:)
    logical,allocatable :: do_pda(:)
    real(dp) :: maxdiff
    real(dp) :: mean_n_photons
    real(dp) :: tolerance

    integer :: ic
    integer :: ipda
    type(grid_cell), allocatable :: pda_cells(:)
    integer,allocatable :: id_pda_cell(:)

    character(len=32) :: envval
    integer :: envlen, envstat
    integer(8) :: clock0, clock1, clock_rate
    logical :: banded_ok, cmp_ok
    real(dp),allocatable :: se_save(:,:), se_band(:,:), em_save(:)
    real(dp),allocatable :: se_stock(:,:), em_band(:), em_stock(:)
    real(dp) :: rel, relmax, res_band, res_stock
    integer :: n_1e3, n_1e2, id

    if(.not.pda_env_read) then
       call get_environment_variable('HYPERION_PDA_SOLVER', envval, envlen, envstat)
       pda_timed = envstat == 0
       pda_banded = pda_timed .and. trim(envval) == 'banded'
       pda_compare = pda_timed .and. trim(envval) == 'compare'
       if(pda_timed) write(*,'(" [pda] HYPERION_PDA_SOLVER = ",A)') trim(envval)
       pda_env_read = .true.
    end if
    call system_clock(clock0, clock_rate)

    mean_n_photons = sum(n_photons) / size(n_photons)

    allocate(specific_energy_prev(geo%n_cells, n_dust))
    allocate(do_pda(geo%n_cells))

    do_pda = n_photons < max(30,ceiling(threshold_pda*mean_n_photons)) .and. sum(density, dim=2) > 0._dp

    call check_allowed_pda(do_pda)

    if(.not.any(do_pda)) then
       write(*,'(" [pda] not necessary for this iteration")')
       return
    end if

    banded_ok = .false.
    if(pda_banded) then
       write(*,'(" [pda] ",I0," PDA cells - using banded direct solver")') count(do_pda)
       tolerance = tolerance_exact
       banded_ok = .true.
       pda_report = .true.
    else if(count(do_pda) < 10000) then
       write(*,'(" [pda] fewer than 10,000 PDA cells - using Gauss pivot method")')
       tolerance = tolerance_exact
    else
       write(*,'(" [pda] more than 10,000 PDA cells - using iterative method")')
       tolerance = tolerance_iter
    end if

    allocate(e_mean(geo%n_cells))

    do ic=1,geo%n_cells
       call update_e_mean(ic)
    end do

    ! Precompute the cells where the PDA will be computed
    allocate(pda_cells(count(do_pda)))

    allocate(id_pda_cell(geo%n_cells))
    id_pda_cell = -1

    ipda = 0
    do ic=1,geo%n_cells
       if(do_pda(ic)) then
          ipda = ipda + 1
          pda_cells(ipda) = new_grid_cell(ic, geo)
          id_pda_cell(ic) = ipda
       end if
    end do

    specific_energy_prev = specific_energy

    do

       specific_energy_prev = specific_energy

       if(pda_compare) then
          if(.not.allocated(se_save)) then
             allocate(se_save(geo%n_cells, n_dust), se_band(geo%n_cells, n_dust), em_save(geo%n_cells))
             allocate(se_stock(geo%n_cells, n_dust), em_band(geo%n_cells), em_stock(geo%n_cells))
             pda_report = .true.
          end if
          se_save = specific_energy
          em_save = e_mean
          call solve_pda_indiv_banded(pda_cells, id_pda_cell, cmp_ok)
          se_band = specific_energy
          em_band = e_mean
          specific_energy = se_save
          e_mean = em_save
       end if

       if(banded_ok) then
          call solve_pda_indiv_banded(pda_cells, id_pda_cell, banded_ok)
          if(.not.banded_ok) then
             ! band matrix too large: stock solver for this solve
             if(count(do_pda) < 10000) then
                tolerance = tolerance_exact
             else
                tolerance = tolerance_iter
             end if
             cycle
          end if
       else if(count(do_pda) < 10000) then
          call solve_pda_indiv_exact(pda_cells, id_pda_cell)
       else
          call solve_pda_indiv_iterative(pda_cells)
       end if

       if(pda_compare .and. cmp_ok) then
          relmax = 0._dp
          n_1e3 = 0
          n_1e2 = 0
          do ipda=1,size(pda_cells)
             ic = pda_cells(ipda)%ic
             do id=1,n_dust
                rel = abs(se_band(ic, id) - specific_energy(ic, id)) / specific_energy(ic, id)
                relmax = max(relmax, rel)
                if(rel > 1.e-3_dp) n_1e3 = n_1e3 + 1
                if(rel > 1.e-2_dp) n_1e2 = n_1e2 + 1
             end do
          end do
          write(*,'(" [pda] banded vs stock: max rel diff ", ES9.2, ", cells > 1e-3: ", I0, ", > 1e-2: ", I0, " of ", I0)') &
               & relmax, n_1e3, n_1e2, size(pda_cells)
          ! residuals of both solutions in the linear system both solvers
          ! set up (coefficients from the specific energy before the solve)
          se_stock = specific_energy
          em_stock = e_mean
          specific_energy = se_save
          e_mean = em_band
          res_band = pda_residual(pda_cells)
          e_mean = em_stock
          res_stock = pda_residual(pda_cells)
          specific_energy = se_stock
          write(*,'(" [pda] max rel residual: banded ", ES9.2, ", stock ", ES9.2)') res_band, res_stock
       end if

       maxdiff = maxval(abs(specific_energy - specific_energy_prev) / specific_energy_prev)

       write(*,'(" [pda] maximum energy difference: ", ES9.2)') maxdiff

       if(maxdiff < tolerance) exit

    end do

    write(*,'(" [pda] converged")')

    if(pda_timed) then
       call system_clock(clock1)
       write(*,'(" [pda] wall time: ", F10.2, " s")') real(clock1 - clock0, dp) / real(clock_rate, dp)
    end if

    deallocate(do_pda)
    deallocate(specific_energy_prev)
    deallocate(e_mean)

    call update_energy_abs_tot()

    call check_energy_abs()

  end subroutine solve_pda

  real(dp) function dtau_rosseland(cell, idir)
    implicit none
    type(grid_cell), intent(in) :: cell
    integer,intent(in) :: idir
    integer :: id
    dtau_rosseland = 0._dp
    do id=1,n_dust
       dtau_rosseland = dtau_rosseland + density(cell%ic,id) * chi_rosseland(id, specific_energy(cell%ic,id)) * cell_width(cell,idir)
    end do
  end function dtau_rosseland

  subroutine solve_pda_indiv_exact(pda_cells, id_pda_cell)

    implicit none

    type(grid_cell),intent(in) :: pda_cells(:)
    integer,intent(in) :: id_pda_cell(:)
    ! Which cells should be used for the PDA

    integer :: direction, wall
    type(grid_cell) :: curr, next
    real(dp) :: dtau_ross_curr, dtau_ross_next, dtau_sum

    integer :: ic
    real(dp) :: coefficient
    real(dp),allocatable :: a(:,:), b(:)
    integer :: id_curr, id_next

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       call update_e_mean(ic)
    end do

    allocate(a(size(pda_cells), size(pda_cells)), b(size(pda_cells)))
    a = 0._dp
    b = 0._dp

    do id_curr=1,size(pda_cells)

       curr = pda_cells(id_curr)

       do wall = 1, geo%n_dim * 2

          direction = int((wall+1)/2)

          next = next_cell(curr, wall)

          dtau_ross_curr = dtau_rosseland(curr, direction)
          dtau_ross_next = dtau_rosseland(next, direction)

          dtau_sum = dtau_ross_curr + dtau_ross_next

          ! If the optical depth is too small, we have to reset it to avoid
          ! issues.
          if(dtau_sum < 1e-100_dp) dtau_sum = 1e-100_dp

          coefficient = 1. / dtau_sum / cell_width(curr, direction)
          coefficient = coefficient * geometrical_factor(wall, curr)

          a(id_curr, id_curr) = a(id_curr, id_curr) - coefficient

          if(id_pda_cell(next%ic) > 0) then
             id_next = id_pda_cell(next%ic)
             a(id_next, id_curr) = coefficient
          else
             b(id_curr) = b(id_curr) - coefficient * e_mean(next%ic)
          end if

       end do

    end do

    call lineq_gausselim(a, b)

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       e_mean(ic) = b(id_curr)
       call update_specific_energy(ic)
    end do

    deallocate(a, b)

  end subroutine solve_pda_indiv_exact

  subroutine solve_pda_indiv_iterative(pda_cells)

    implicit none

    type(grid_cell),intent(in) :: pda_cells(:)
    ! Which cells should be used for the PDA

    integer :: direction, wall
    type(grid_cell) :: curr, next
    real(dp) :: dtau_ross_curr, dtau_ross_next

    real(dp) :: coefficient
    real(dp) :: a, b
    integer :: id_curr, ic

    real(dp) :: max_e_diff, e_diff, e_new

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       call update_e_mean(ic)
    end do

    do

       max_e_diff = 0.
       do id_curr=1,size(pda_cells)

          curr = pda_cells(id_curr)

          a = 0._dp
          b = 0._dp

          do wall = 1, geo%n_dim * 2

             direction = int((wall+1)/2)

             next = next_cell(curr, wall)

             dtau_ross_curr = dtau_rosseland(curr, direction)
             dtau_ross_next = dtau_rosseland(next, direction)

             coefficient = 1. / (dtau_ross_curr + dtau_ross_next) / cell_width(curr, direction)
             coefficient = coefficient * geometrical_factor(wall, curr)

             a = a - coefficient
             b = b - coefficient * e_mean(next%ic)

          end do

          e_new = b/a

          e_diff = abs(e_new - e_mean(curr%ic)) / e_mean(curr%ic)
          if(e_diff > max_e_diff) max_e_diff = e_diff

          e_mean(curr%ic) = e_new

       end do

       if(max_e_diff < tolerance_iter) exit

    end do

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       call update_specific_energy(ic)
    end do

  end subroutine solve_pda_indiv_iterative

  subroutine solve_pda_indiv_banded(pda_cells, id_pda_cell, ok)

    ! The linear system of solve_pda_indiv_exact (one equation per PDA cell:
    ! the sum over its walls of coefficient * (e_next - e_curr) is zero, with
    ! the terms of non-PDA neighbours on the right-hand side; a wall that
    ! wraps onto the cell itself contributes nothing), solved directly as a
    ! band matrix.  The cells are numbered with the first or the second grid
    ! index running fastest, whichever gives the smaller bandwidth, and the
    ! band LU factorisation runs without pivoting: every diagonal is minus the
    ! sum of the cell's wall coefficients and each off-diagonal is one of
    ! them, so the matrix is diagonally dominant, strictly so in the rows of
    ! cells next to non-PDA cells, which every connected PDA region has
    ! (check_allowed_pda excludes the grid edges).  Cost ~ n * bandwidth**2
    ! instead of ~ n**3 (exact) or many Gauss-Seidel sweeps (iterative).
    ! ok is set to .false. without solving if the band matrix would exceed
    ! pda_band_max_elements.

    implicit none

    type(grid_cell),intent(in) :: pda_cells(:)
    integer,intent(in) :: id_pda_cell(:)
    logical,intent(out) :: ok

    integer :: n, nw, i, j, k, kmax, r, bw, bw1, bw2, wall, direction, ic, id_next
    integer :: i1, i2, i3
    integer, allocatable :: nbr(:,:), perm(:), inv(:)
    real(dp), allocatable :: coef(:,:), rhs(:), ab(:,:), x(:)
    type(grid_cell) :: curr, next
    real(dp) :: dtau_ross_curr, dtau_ross_next, dtau_sum, coefficient, f

    n = size(pda_cells)
    nw = geo%n_dim * 2

    do i=1,n
       ic = pda_cells(i)%ic
       call update_e_mean(ic)
    end do

    ! wall coefficients (as in solve_pda_indiv_exact) and PDA neighbours
    ! (index into pda_cells, 0 for a non-PDA cell)
    allocate(nbr(nw, n), coef(nw, n), rhs(n))
    rhs = 0._dp
    do i=1,n
       curr = pda_cells(i)
       do wall=1,nw
          direction = int((wall+1)/2)
          next = next_cell(curr, wall)
          dtau_ross_curr = dtau_rosseland(curr, direction)
          dtau_ross_next = dtau_rosseland(next, direction)
          dtau_sum = dtau_ross_curr + dtau_ross_next
          if(dtau_sum < 1e-100_dp) dtau_sum = 1e-100_dp
          coefficient = 1. / dtau_sum / cell_width(curr, direction)
          coefficient = coefficient * geometrical_factor(wall, curr)
          coef(wall, i) = coefficient
          id_next = id_pda_cell(next%ic)
          if(id_next > 0) then
             nbr(wall, i) = id_next
          else
             nbr(wall, i) = 0
             rhs(i) = rhs(i) - coefficient * e_mean(next%ic)
          end if
       end do
    end do

    ! numbering 1: pda_cells order (cell id order, first index fastest)
    bw1 = 0
    do i=1,n
       do wall=1,nw
          if(nbr(wall, i) > 0) bw1 = max(bw1, abs(nbr(wall, i) - i))
       end do
    end do

    ! numbering 2: second index fastest; perm(new) = old, inv(old) = new
    allocate(perm(n), inv(n))
    k = 0
    do i3=1,geo%n3
       do i1=1,geo%n1
          do i2=1,geo%n2
             ic = cell_id(i1, i2, i3, geo)
             if(id_pda_cell(ic) > 0) then
                k = k + 1
                perm(k) = id_pda_cell(ic)
             end if
          end do
       end do
    end do
    if(k /= n) call error("solve_pda_indiv_banded", "PDA cell count mismatch")
    do k=1,n
       inv(perm(k)) = k
    end do
    bw2 = 0
    do i=1,n
       do wall=1,nw
          if(nbr(wall, i) > 0) bw2 = max(bw2, abs(inv(nbr(wall, i)) - inv(i)))
       end do
    end do

    if(bw1 <= bw2) then
       do i=1,n
          inv(i) = i
       end do
       bw = bw1
    else
       bw = bw2
    end if

    if(real(2 * bw + 1, dp) * real(n, dp) > pda_band_max_elements) then
       write(*,'(" [pda] band matrix ",I0," x ",I0," too large - using the stock solver")') n, 2 * bw + 1
       ok = .false.
       return
    end if
    ok = .true.

    if(pda_report) then
       write(*,'(" [pda] banded solver: bandwidth ",I0," (first index fastest ",I0,", second ",I0,")")') bw, bw1, bw2
       pda_report = .false.
    end if

    ! band storage: ab(j - r, r) = A(r, j) for equation r, unknown j (new numbering)
    allocate(ab(-bw:bw, n), x(n))
    ab = 0._dp
    do i=1,n
       r = inv(i)
       x(r) = rhs(i)
       do wall=1,nw
          j = nbr(wall, i)
          if(j == i) cycle
          ab(0, r) = ab(0, r) - coef(wall, i)
          if(j > 0) ab(inv(j) - r, r) = ab(inv(j) - r, r) + coef(wall, i)
       end do
    end do

    ! LU factorisation without pivoting, applied to the right-hand side
    do k=1,n-1
       kmax = min(n, k + bw)
       do r=k+1,kmax
          if(ab(k - r, r) /= 0._dp) then
             f = ab(k - r, r) / ab(0, k)
             ab(k + 1 - r:kmax - r, r) = ab(k + 1 - r:kmax - r, r) - f * ab(1:kmax - k, k)
             x(r) = x(r) - f * x(k)
             ab(k - r, r) = 0._dp
          end if
       end do
    end do

    ! back substitution
    do r=n,1,-1
       kmax = min(n, r + bw)
       if(kmax > r) x(r) = x(r) - sum(ab(1:kmax - r, r) * x(r + 1:kmax))
       x(r) = x(r) / ab(0, r)
    end do

    do i=1,n
       ic = pda_cells(i)%ic
       e_mean(ic) = x(inv(i))
       call update_specific_energy(ic)
    end do

    deallocate(nbr, coef, rhs, perm, inv, ab, x)

  end subroutine solve_pda_indiv_banded

  real(dp) function pda_residual(pda_cells)

    ! Largest relative residual |sum_w c_w (e_next - e_curr)| / sum_w c_w e_next
    ! over the PDA cells, for the current e_mean, with the wall coefficients of
    ! solve_pda_indiv_exact (compare mode only).

    implicit none

    type(grid_cell),intent(in) :: pda_cells(:)
    integer :: id_curr, wall, direction
    type(grid_cell) :: curr, next
    real(dp) :: dtau_sum, coefficient, a, b

    pda_residual = 0._dp
    do id_curr=1,size(pda_cells)
       curr = pda_cells(id_curr)
       a = 0._dp
       b = 0._dp
       do wall = 1, geo%n_dim * 2
          direction = int((wall+1)/2)
          next = next_cell(curr, wall)
          if(next%ic == curr%ic) cycle
          dtau_sum = dtau_rosseland(curr, direction) + dtau_rosseland(next, direction)
          if(dtau_sum < 1e-100_dp) dtau_sum = 1e-100_dp
          coefficient = 1. / dtau_sum / cell_width(curr, direction) * geometrical_factor(wall, curr)
          a = a + coefficient
          b = b + coefficient * e_mean(next%ic)
       end do
       pda_residual = max(pda_residual, abs(b - a * e_mean(curr%ic)) / abs(b))
    end do

  end function pda_residual

end module grid_pda
