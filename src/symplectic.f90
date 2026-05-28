module symplectic
  implicit none
  private
  public :: is_symplectic, symplectic_compose, random_symplectic

contains

  logical function is_symplectic(mat, nn, tol)
    real(8), intent(in) :: mat(:,:)
    integer, intent(in) :: nn
    real(8), intent(in) :: tol
    real(8), allocatable :: jj(:,:), work(:,:), check(:,:)
    integer :: sz, ii, jjj

    sz = 2 * nn
    allocate(jj(sz,sz), work(sz,sz), check(sz,sz))

    ! Build canonical J = [0 I; -I 0]
    jj = 0.0d0
    do ii = 1, nn
      jj(ii, nn + ii) = 1.0d0
      jj(nn + ii, ii) = -1.0d0
    end do

    work = matmul(transpose(mat), jj)
    check = matmul(work, mat)

    is_symplectic = .true.
    do jjj = 1, sz
      do ii = 1, sz
        if (abs(check(ii, jjj) - jj(ii, jjj)) > tol) then
          is_symplectic = .false.
          return
        end if
      end do
    end do

    deallocate(jj, work, check)
  end function

  subroutine symplectic_compose(mat1, mat2, res, nn)
    real(8), intent(in) :: mat1(:,:), mat2(:,:)
    real(8), intent(out) :: res(:,:)
    integer, intent(in) :: nn

    res = matmul(mat1, mat2)
  end subroutine

  subroutine random_symplectic(mat, nn, seedval)
    real(8), intent(out), allocatable :: mat(:,:)
    integer, intent(in) :: nn, seedval
    real(8), allocatable :: ss(:,:), js(:,:), term(:,:)
    integer :: sz, ii, jj, kk, cnt
    real(8) :: harvest
    integer :: seed_arr(8)

    sz = 2 * nn
    allocate(mat(sz,sz), ss(sz,sz), js(sz,sz), term(sz,sz))

    ! Seed RNG
    call random_seed(size=cnt)
    if (cnt > 8) cnt = 8
    do ii = 1, cnt
      seed_arr(ii) = seedval + ii
    end do
    call random_seed(put=seed_arr)

    ! Build canonical J into mat temporarily
    mat = 0.0d0
    do ii = 1, nn
      mat(ii, nn + ii) = 1.0d0
      mat(nn + ii, ii) = -1.0d0
    end do

    ! Random symmetric matrix
    ss = 0.0d0
    do ii = 1, sz
      do jj = ii, sz
        call random_number(harvest)
        ss(ii, jj) = (harvest - 0.5d0) * 0.5d0
        ss(jj, ii) = ss(ii, jj)
      end do
    end do

    js = matmul(mat, ss)

    ! Taylor exp(JS), 12 terms
    mat = 0.0d0
    term = 0.0d0
    do ii = 1, sz
      mat(ii, ii) = 1.0d0
      term(ii, ii) = 1.0d0
    end do

    do kk = 1, 12
      term = matmul(term, js) / dble(kk)
      mat = mat + term
    end do

    deallocate(ss, js, term)
  end subroutine

end module
