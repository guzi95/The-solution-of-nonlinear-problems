! Operacoes de Algebra Linear. Ferramentas necessarias para o processo de inversao.
! Autora: JULIANA DA SILVA FURTADO - Brasil, Para

MODULE levenberg_marq_module
IMPLICIT NONE
save
CONTAINS

!----------------------------------------------------------------------------------------
! Calcula os valores de gravidade a partir dos parâmetros x e o erro entre calculado e observado C.
subroutine gravity_and_c_values(n,t,x,y,obs,g,c)
implicit none
real, parameter:: pi= 3.141, k= 6.672, d= 0.25
integer, intent(in):: n, t
real, intent(in):: x(t), y(n), obs(n)
real, intent(out):: g(n), c(n)

    g(:)= (((4*pi)/3)*(k*d*x(1)*(x(2)**3)))/(((y(:)**2)+(x(1)**2))**(3/2))
    c(:)= obs(:)-g(:)

end subroutine

!------------------------------------------------------------------------------------
! Calcula as matrizes atc para cada novo x.
subroutine matrizes_ata_atc (n,t,lambda,x,y,c,a,atc,ata) 
implicit none
real, parameter:: pi= 3.141, k= 6.672, d= 0.25
real, intent(in):: lambda, x(t), y(n), c(n)
real, intent(out):: a(n,t), atc(t), ata(t,t)
integer, intent(in):: n, t    
integer:: i, j

    !! Derivada parcial em função da variavel "profundidade".
    a(:,1)= ((4*pi)*k*d*(x(2)**3)*((y(:)**2)-(2*(x(1)**2))))/(((y(:)**2)+(x(1)**2))**(5/2))
    !! Derivada Parcial em função da variavel "raio".
    a(:,2)= (4*pi*k*d*(x(2)**2)*x(1))/(((y(:)**2)+(x(1)**2))**(3/2))
    
    do i= 1, t
        do j= 1, t
            ata(i,j)= dot_product(a(:,i),a(:,j))
        end do
        ata(i,i)= ata(i,i)+lambda
        atc(i)= dot_product(a(:,i),c(:))
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
subroutine sub2_solve (n,a,L,atc,y)  
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

! Ultimo elemento de x.
y(n)= y(n) / A(n,n)

do i= n-1, 1, -1
plus= 0.0
    do j= i+1, n
    plus= (A(i,j)*y(j)) + plus
    enddo
    y(i)= (y(i) - plus) / A(i,i)
enddo 
 
end subroutine sub2_solve

end module
