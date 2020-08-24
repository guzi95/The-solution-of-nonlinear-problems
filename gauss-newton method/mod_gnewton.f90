MODULE gauss_newton_method
implicit none
save
CONTAINS

!--------------------------------------------------------------------------------------------
! Calcula os valores de g para cada novo x.
subroutine gravity_values(n,t,g,x,y,o,c)
implicit none
integer, intent(in):: n, t
real, parameter:: pi= 3.141593, k= 6.672, d= 0.25
real, intent(in):: x(t), y(n), o(n)
real, intent(out):: g(n),c(n)

    g(:)= ((((4*pi)/3)*k*d*(x(1)**3)*x(2))/((y(:)**2)+((x(2)**2)**(3/2))))
    c(:)= o(:)-g(:)
    
end subroutine

!---------------------------------------------------------------------------------------------
! Calcula as matrizes A, At, At.A e At.c para cada novo x.
subroutine at_matrix (n,t,a,atc,ata,c,x,y) 
implicit none
real, parameter:: pi= 3.141593, k= 6.672, d= 0.25
integer, intent(in):: n, t
real, intent(in):: x(t),y(n),c(n)
real, intent(out):: a(n,t), ata(t,t), atc(t)
integer:: i, j

! Calcula os elementos de A.
    A(:,1)= (4*pi*k*d*(x(1)**2)*x(2))/(((y(:)**2)+(x(2)**2))**(3/2))
    A(:,2)= (((4*pi)*k*d*(x(1)**3)*((y(:)**2)-(2*(x(2)**2))))/(((y(:)**2)+(x(2)**2))**(5/2)))

! Calcula AtA e At.c.
do i= 1, t
    do j= 1, t
        ata(i,j)= dot_product(a(:,i),a(:,j))
    end do
    atc(i)= dot_product(a(:,i),c(:)) + (2/(n-2))
end do

end subroutine

!---------------------------------------------------------------------------------------------
! Decompoe a matriz At.A
subroutine sub1_decomp (n,a,l)
implicit none
integer, intent(in):: n
real, intent(inout):: a(n,n)
real, intent(out):: l(n,n)
integer:: i, p
    
do p= 1, n-1
    do i= p+1, n
        a(i,p)= a(i,p)/a(p,p)
        a(i,p+1:n)= a(i,p+1:n) - (l(i,p)*a(p,p+1:n))
    enddo
enddo

end subroutine sub1_decomp

!---------------------------------------------------------------------------------------------
! Resolve o problema inverso AtA.(delta)= Atc
subroutine sub2_solve (n,A,L,atc,y)  
implicit none
integer, intent(in):: n
real, intent(in):: a(n,n), l(n,n), atc(n)
real, intent(out):: y(n)
integer:: i, j
real::plus
    
y(1)= atc(1)  
do i= 2, n
    plus= 0.0
    do j= 1, n-1
        plus= plus + (L(i,j)*y(j))
    enddo
    y(i)= atc(i) - plus
enddo

! ultimo elemento de x
y(n)= y(n) / A(n,n)

do i= n-1, 1, -1
plus= 0.0
    do j= i+1, n
    plus= (A(i,j)*y(j)) + plus
    enddo
    y(i)= (y(i) - plus) / A(i,i)
enddo 
 
end subroutine sub2_solve

   
END MODULE gauss_newton_method
