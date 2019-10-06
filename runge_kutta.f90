module runge_kutta
    implicit none
    public order_4_step
    interface
        real(8) function ode_interface(t, y, n)
        integer, intent(in) :: n
        real(8), intent(in) :: t
        real(8), intent(in) :: y(n)
        end function ode_interface
    end interface
contains
function order_4_step(f, t, y, h, n) result(y_new)
    procedure(ode_interface) :: f
    integer, intent(in) :: n
    real(8), intent(in) :: t
    real(8), intent(in) :: h
    real(8), intent(in) :: y(n)
    
    real(8) :: k(4,n)
    real(8) :: y_new(n)
    
    k(1,:) = h/6.d0*f(t, y, n)
    k(2,:) = h/3.d0*f(t + 1.d0/2.d0*h, y + 1.d0/2.d0*k(1,:), n)
    k(3,:) = h/3.d0*f(t + 1.d0/2.d0*h, y + 1.d0/2.d0*k(2,:), n)
    k(4,:) = h/6.d0*f(t + h, y + k(3,:), n)
    y_new(1:n) = y(1:n) + sum(k(:,1:n))
    return
end function order_4_step
end module runge_kutta