! Operacoes de Algebra Linear. Ferramentas necessarias para o processo de inversao.
! Autora: JULIANA DA SILVA FURTADO - Brasil, Para

MODULE steepest_descent_method
implicit none
save
CONTAINS

! Calcula os valores de gravidade a partir dos parâmetros x e o erro entre dados calculados e observados C.
subroutine gravity_and_c_values(n,t,x,y,obs,g,c)
implicit none
real, parameter:: pi= 3.141, k= 6.672, d= 0.25
integer, intent(in):: n, t
real, intent(in):: x(t), y(n), obs(n)
real, intent(out):: g(n), c(n)

    g(:)= (((4*pi)/3)*(k*d*x(1)*(x(2)**3)))/(((y(:)**2)+(x(1)**2))**(3/2))
    c(:)= obs(:)-g(:)

end subroutine

! Calcula as matrizes atc para cada novo x.
subroutine atc_matrix(n,t,u,x,y,c,grad,delta,a,atc) 
implicit none
real, parameter:: pi= 3.141, k= 6.672, d= 0.25
real:: u
real, intent(in):: x(t), y(n), c(n)
real, intent(out):: a(n,t), atc(t), grad(t), delta(t)
integer, intent(in):: n, t    
integer:: i

    !! Derivada parcial em função da variável "profundidade".
    A(:,1)= ((4*pi)*k*d*(x(2)**3)*((y(:)**2)-(2*(x(1)**2))))/(((y(:)**2)+(x(1)**2))**(5/2))
    !! Derivada Parcial em função da variável "raio".
    A(:,2)= (4*pi*k*d*(x(2)**2)*x(1))/(((y(:)**2)+(x(1)**2))**(3/2))

do i= 1, t
    atc(i)= dot_product(a(:,i),c(:))
end do

    ! Gradiente.
    grad(:)= atc(:)
    ! "passo"
    delta(:)= grad(:)/(u)
    
end subroutine

end module
