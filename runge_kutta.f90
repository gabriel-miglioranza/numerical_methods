module runge_kutta
    implicit none
    public order_4_step, ode_interface
    interface
        function ode_interface(t, y, rpar, ipar)
            integer, intent(in), optional :: ipar(:)
            real(8), intent(in) :: t
            real(8), intent(in) :: y(:)
            real(8), intent(in), optional :: rpar(:)
            real(8) :: ode_interface(size(y)) 
        end function ode_interface
    end interface
contains
function order_4_step(f, t, y, h, rpar, ipar)
    procedure(ode_interface) :: f
    integer, intent(in), optional :: ipar(:)
    real(8), intent(in) :: t
    real(8), intent(in) :: h
    real(8), intent(in) :: y(:)
    real(8), intent(in), optional :: rpar(:)
    real(8) :: order_4_step(size(y)) 
    real(8) :: k(4,size(y))
    integer :: i
    k(1,:) = h/6.d0*f(t, y, rpar, ipar)
    k(2,:) = h/3.d0*f(t + 1.d0/2.d0*h, y + 1.d0/2.d0*k(1,:), rpar, ipar)
    k(3,:) = h/3.d0*f(t + 1.d0/2.d0*h, y + 1.d0/2.d0*k(2,:), rpar, ipar)
    k(4,:) = h/6.d0*f(t + h, y + k(3,:), rpar, ipar)
    do i = 1, size(y)
        order_4_step(i) = y(i) + sum(k(:,i))
    end do
    return
end function order_4_step
end module runge_kutta