
program main_preg1

use nrtype
use nr
use utils
use general_fit

implicit none

integer, parameter                     :: ncol=3, npts=100
integer                                :: n, i
real(8),dimension(npts)                :: theta
real(8),dimension(:,:),allocatable     :: data
real(8), dimension(0:5)                :: a1, a2
real(8), dimension(0:5)                :: sigma_a1, sigma_a2
real(8)                                :: q1, q2,theta_min,theta_max

! Leyendo datos de 'gamma.dat'

call readncol(ncol,data,n)

theta_min = pi_d/2.0_dp
theta_max = pi_d

theta = linspace(theta_min,theta_max,npts)

! base legendre
call linear_fit(fbase1,data(:,1),data(:,2),data(:,3),a1,sigma_a1,q1)

write(*,'(3X,A)') 'Coeficientes con el ajuste de polinomios de Legendre:'
print*
write(*,'(2F15.7)') (a1(i),sigma_a1(i) , i=0,size(a1)-1)
print*
write(*,*) 'q = ',q1

open(100,file='base_legendre.dat')
write(100,'(2e15.7)') (cos(theta(i)), Dot_product(fbase1(cos(theta(i))),a1), i=1,npts)
close(100)

print*
print*, "Plotting results with a gnuplot script ..."
call system("gnuplot base_legendre.plt")
read*

!base cosenos  ##########

call linear_fit(fbase2,data(:,1),data(:,2),data(:,3),a2,sigma_a2,q2)

write(*,'(3X,A)') 'Coeficientes con el ajuste de base de cosenos:'
print*
write(*,'(2F15.7)') (a2(i),sigma_a2(i) , i=0,size(a2)-1)
print*
write(*,*) 'q = ',q2

open(200,file='base_cos.dat')
write(200,'(2e15.7)') (cos(theta(i)), Dot_product(fbase2(cos(theta(i))),a2), i=1,npts)
close(200)

print*
print*, "Plotting results with a gnuplot script ..."
call system("gnuplot base_cos.plt")
read*



contains

    function fbase1(x)
    real(8),intent(in)    :: x
    real(8),dimension(0:5)  :: fbase1
    real(8),dimension(0:10) :: legen_coef
    integer                 :: n,i

    legen_coef(0) = 1
    legen_coef(1) = x

    do n=2,size(legen_coef)-1
        legen_coef(n) =(1_8/real(n,8))*(real(2*n-1,8)*x*legen_coef(n-1) - real(n-1,8)*legen_coef(n-2))
    end do

    do i=0,size(fbase1)-1
        fbase1(i) = legen_coef(2*i)
    end do

    end function
!#######################################################
    function fbase2(x)
    real(8),intent(in)    :: x
    real(8),dimension(0:5)  :: fbase2
    integer                 :: i

    do i=0,size(fbase2)-1
        fbase2(i) = x**(2*i)
    end do

    end function

end program