module gauss_seidel
    implicit none
    public solve
    private gauss_seidel_iter
contains
function solve(a, x_initial, b, n, atol, rtol, relaxation) result(x_new)
    integer, intent(in) :: n
    real(8), intent(in) :: atol, rtol
    real(8), intent(in) :: x_initial(n)
    real(8), intent(in) :: b(n)
    real(8), intent(in) :: a(n, n)
    real(8), intent(in), optional :: relaxation

    integer :: i
    real(8) ::  x_new(n), x_old(n), w = 1.d0
    
    if (present(relaxation)) w = relaxation
    x_new = x_initial
    do
        do i = 1, n
        x_old = x_new
        x_new(i) = x_old(i) + w*gauss_seidel_iter(a, x_new, x_old, b, n, i) 
        end do
        if(maxval(abs(x_new-x_old)) < atol + rtol * maxval(abs(x_new))) exit
    end do
    return
end function solve
pure real(8) function gauss_seidel_iter(a, x_new, x_old, b, n, i) 
    integer, intent(in) :: n
    integer, intent(in) :: i
    real(8), intent(in) :: x_new(n)
    real(8), intent(in) :: x_old(n)
    real(8), intent(in) :: b(n)
    real(8), intent(in) :: a(n, n)

    gauss_seidel_iter = 1.d0/a(i,i) * (b(i) - & 
                    sum(a(i,:i-1)*x_new(:i-1)) &
                    - sum(a(i,i+1:)*x_old(i+1:)))
    return
end function gauss_seidel_iter
end module gauss_seidel