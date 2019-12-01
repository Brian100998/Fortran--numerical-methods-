
module general_fit
use nr
use nrtype
use nrutil
use utils
implicit none

contains
!#######################################################################
    function desing_matrix(fbase,xi,sigma_i)

    real(8), intent(in), dimension(:)    :: xi
    real(8), intent(in), dimension(:)    :: sigma_i
    real(8), dimension(:,:), allocatable :: desing_matrix 
    integer                              :: i,j,n,m 
    interface    
        function fbase(x)
        implicit none
        real(8),intent(in)     :: x
        real(8),dimension(0:5) :: fbase
        end function
    end interface


    n = size(xi)
    m = 6

    allocate(desing_matrix(n,m))

    do i=1,n
        desing_matrix(i,:) = fbase(xi(i))/sigma_i(i)
    end do

    end function
!#######################################################################
    subroutine linear_fit(fbase,xi,yi,sigma_i,a,sigma_a,q)


    real(8), intent(in), dimension(:)    :: xi
    real(8), intent(in), dimension(:)    :: yi
    real(8), intent(in), dimension(:)    :: sigma_i
    real(8), intent(out), dimension(0:5,1) :: a
    real(8), intent(out), dimension(0:5) :: sigma_a
    real(8), intent(out)                 :: q
    real(8), dimension(:,:), allocatable :: matrix_coeff,A_desig,b,vector_b,chi,a_pri,chi_prim
    real(8)                              :: chi2,sum_int
    integer                              :: y_siz,p,m,n,i,k,j,l,h
    interface    
        function fbase(x)
        implicit none
        real(8),intent(in)     :: x
        real(8),dimension(0:5) :: fbase
        end function
    end interface   
    ! matrix b
    y_siz = size(yi)

    allocate(b(y_siz,1))

    do p=1,y_siz
        b(p,1) = yi(p)/sigma_i(p)
    end do
    
    ! a coefficients

    allocate(A_desig(size(xi),6))


    A_desig = desing_matrix(fbase,xi,sigma_i)

    m = size(A_desig(1,:))
    n = size(A_desig(:,1))


    allocate(matrix_coeff(m,m))
    allocate(vector_b(m,1))
    matrix_coeff = Matmul(Transpose(A_desig), A_desig)
    vector_b = Matmul(Transpose(A_desig),b)



    call gaussj(matrix_coeff,vector_b)

    a = vector_b

    ! sigma_a elements, guassj devolvio  matrix_coeff = (A^t * A)^-1

    do i=0,size(a)-1
        sigma_a(i) = sqrt(matrix_coeff(i+1,i+1))
    end do

    ! falta determinar q

    allocate(chi(n,1))
    allocate(a_pri(1,size(a(:,1))))
    a_pri = Transpose(a)

    do j=1,n
        chi(j,1) = b(j,1) - dot_product(a_pri(1,:),A_desig(j,:))
    end do

    allocate(chi_prim(1,n))
    chi_prim = Transpose(chi)
    chi2 = dot_product(chi_prim(1,:),chi_prim(1,:))

    if(N>2) q = gammq(0.5_dp*real(n-m,8) ,0.5_dp*chi2)
    !q=1
    end subroutine
!####################################################################

end module
