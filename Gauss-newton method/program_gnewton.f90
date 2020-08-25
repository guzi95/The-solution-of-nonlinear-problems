! Este programa calcula um problema inverso pelo metodo de Gauss-Newton.
! O problema: Em posse dos dados observado em campo (aqui sao componentes gravitacionais), deseja-se estimar a
! profundidade e o raio da anomalia geologica esferica em subsuperficie atraves do metodo de Gauss-Newton.
! Autora: JULIANA DA SILVA FURTADO - Brasil, Para

program GaussNewton_Method
use gauss_newton_method
implicit none
real:: s
real, allocatable:: y(:), x(:), o(:), a(:,:), at(:,:), atc(:), ata(:,:), c(:), g(:), L(:,:), delta(:)
integer:: i, j, n, t        ! Numero de estacoes e parametros.

! Quantidade de estacoes e de parametros.
open (file='estacoes_e_parametros.dat', unit=1, status='old', action='read')
read(1,*) n, t
close(1)

allocate (y(n), x(t), o(n), a(n,t), at(t,n), atc(t), ata(t,t), c(n), g(n), L(t,t), delta(t))

! Posicao das estacoes.
open (unit= 2, file='estacoes.dat', status='old', action='read')
read(2,*) y(:)
close(2)

! Dados obervados em campo.
open (unit= 3, file='observado3.dat', status='old', action='read')
read(3,*) o(:)
close(3)

! Raio e profundidade iniciais ("chute")
open (unit= 4, file='chute_x3.dat', status='old', action='read')
read(4,*) x(:)
close(4)

! "Passo" em direção ao minimo.
open(unit= 5, file= 'delta.dat', status= 'old', action= 'write')

! Valores dos parâmetros(x).
open(unit= 6, file= 'final_x.dat', status= 'old', action= 'write')

! Escreva os valores da função s(merit function).
open(unit= 7, file= 'merit_function.dat', status= 'old', action= 'write')

! Residuos.
open(unit= 8, file='c_residuo.dat', status='old', action='write')

! Valores de gravidade calculados
open(unit= 9, file='calculado.dat', status='old', action='write')

! Matriz A.
open(unit= 10, file='a.dat', status='old', action='write')

! Matriz AtC.
open(unit= 11, file='atc.dat', status='old', action='write')

! Inicia o processo iterativo para a estimativa dos parametros.
delta(:)= 0.
do
    x(:)= x(:) + delta(:)    ! Novo valor de x.
    write(6,*) (x(i),'km', i= 1, t)
    
    CALL gravity_values(n,t,g,x,y,o,c)
    
    ! Escreve os valores de g
    do i= 1, n
        write(9,*) g(i)
    end do
    
    ! Escreve a matriz C
    do i= 1, n
        write(8,*) c(i)
    end do
    
    ! Funcao "merit". Soma do quadrado dos residuos.
    s= sum(c**2)/(n-2)
    write(7,*) s
    
    CALL at_matrix (n,t,a,atc,ata,c,x,y)  
    
    ! Escreve a matriz A.
    do i= 1, n
        write(10,*) (a(i,j), j= 1, t)
    end do
    
    ! Escreve a matriz AtC.
    write(11,*) (atc(j), j= 1, t)
    
    CALL sub1_decomp (t,ata,L)
    
    CALL sub2_solve(t,ata,L,atc,delta)    ! Novo valor de delta.
    write(5,*) delta(:)
    
    ! Condicao para interromper a iteracao.
    ! Se as condicoes nao forem satisfeitas entao calcule o novo valor de x para uma nova iteracao. 
    if ((abs(delta(1))<=+1.0E-5) .AND. (abs(delta(2))<=+1.0E-5)) exit

end do

close(5)
close(7)
close(8)
close(9)
close(10)
close(11)

write(6,*) 'raio e profundidade, respectivamente:'
do i= 1, t
    write(6,100) x(i)
    100 format(' ',1F10.4)
end do

close(6)
deallocate(y,x,o,a,at,ata,c,g,L,delta)

end program
