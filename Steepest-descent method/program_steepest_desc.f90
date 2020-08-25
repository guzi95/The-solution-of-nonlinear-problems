! Este programa calcula o problema inverso pelo metodo do Maximo Declive. 
! Dos dados observados em campo (neste caso, sao componentes gravitacionais) deseja-se
! estimar a profundidade (parametro 1) e o raio (parametro 2) de uma anomalia geologica esferica
! que esta na subsuperficie terrestre.
! Autora: JULIANA DA SILVA FURTADO - Brasil, Pará

program program_steepest_descent
use steepest_descent_method
implicit none
real, parameter:: tol= 1.0E-5
real:: u= 0., s= 0.
real, allocatable:: y(:), obs(:), initial_x(:), x(:), c(:), grad(:), delta(:), atc(:), g(:), a(:,:)
integer:: n, t, i, j   

! Quantidade de estacoes, parametros e constante u.
open (unit=1, file='estacoes_e_parametros.dat', action='read', status='old')
 read(1,*) n, t
close(1)

! Valor de mi (u).
open (unit= 14, file='u1.dat', action= 'read', status= 'old')
read(14,*) u
close(14)

allocate (y(n),obs(n),initial_x(t),x(t),c(n),grad(t),delta(t),a(n,t),atc(t),g(n))

! Posicoes das estacoes.
open (unit= 2, action='read', file='estacoes.dat', status='old')
read(2,*) y(:)
close(2)

! Dados observados em campo (dados sinteticos).
open (unit= 3, action='read', file='observado1.dat', status='old')
read(3,*) obs(:)
close(3)

! Valores iniciais dos parametros ("chute").
open(unit= 4, action= 'read', file= 'kick_x1.dat', status= 'old')
read(4,*) initial_x(:)
close(4)

! valores finais dos parametros.
open(unit= 5, file= 'final_x.dat', status='old', action='write')

! Residuos.
open(unit= 6, file= 'c_residuo.dat', status='old', action='write')

! gradiente.
open(unit= 7, file='grad.dat', status='old', action='write')

! Delta.
open(unit= 8, file='delta.dat', status='old', action='write')

! valores de gravidade calculados
open(unit= 9, file='calculado.dat', status='old', action='write')

! Matriz A.
open(unit= 10, file='a.dat', status='old', action='write')

! Matriz At.C.
open(unit= 11, file='atc.dat', status='old', action='write')

! Write out s(merit function).
open(unit= 12, file= 'merit_function.dat', status= 'old', action= 'write')

! Inicia as estimativas dos parametros profundidade e raio da anomalia geologica esferica.

do
    ! Calcula o campo gravitacional (g) obtido dos valores do parametros, e o erro (c) entre esses valores e os dados observados.
    CALL gravity_and_c_values(n,t,initial_x,y,obs,g,c)
    
    ! Write out g.
    do i= 1, n
        write(9,*) g(i)
    end do
    
    ! Write out c.
    do i= 1, n
        write(6,*) c(i)
    end do
    
    ! Merit function. The sum of the squares of the residues.
    s= sum(c**2)/(n-2)
    write(12,*) s
    
    ! Calcula o gradiente (grad) and the passo (delta) que sera dado em direcao ao minimo da funcao.
    CALL atc_matrix(n,t,u,initial_x,y,c,grad,delta,a,atc)
    
    ! Write out A matrix.
    do i= 1, n
        write(10,*) (a(i,j), j= 1, t)
    end do
    
    ! Escreva a matriz At.c.
    write(11,*) (atc(j), j= 1, t)
    
    ! Escreva os valores de gradiente.
    write(7,*) (grad(i), i= 1, t)
    
    ! Escreva os valores de delta.
    write(8,*) (delta(i), i= 1, t)
    
    ! x(p+1)
    x(:)= initial_x(:) + delta(:)   
    write(5,*) (x(i), i= 1, t)
    
    ! Condicao para parar o processo iterativo.
    if ((abs(delta(1)) <= +tol) .AND. (abs(delta(2)) <= +tol)) exit
    ! If the condition is not met:
    initial_x(:)= x(:)

end do

! Write out x final values.
write(5,*) '# Depth(z) and radius(a) estimated values, in kilometers (km), respectively:', (x(i), i= 1, t)

close(5)
close(6)
close(7)
close(8)
close(9)
close(10)
close(11)
close(12)
close(13)

deallocate (y,obs,initial_x,x,c,grad,delta,atc)

end program
