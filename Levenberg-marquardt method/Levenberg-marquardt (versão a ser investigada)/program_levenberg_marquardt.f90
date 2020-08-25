! Este programa calcula o problema inverso pelo metodo de Levenberg-Marquardt. 
! Dos dados observados em campo (neste caso, sao componentes gravitacionais) deseja-se
! estimar a profundidade (parametro 1) e o raio (parametro 2) de uma anomalia geologica 
! em formato de esfera que esta na subsuperficie.
! Autor: JULIANA DA SILVA FURTADO - Brasil, Para

program Levenberg_Marquardt_Method
use levenberg_marquardt_module
implicit none
real:: tol= 1.0E-5, k= 2, lambda, s, previous_s, it
real, allocatable:: y(:), obs(:), previous_x(:), x(:), c(:), delta(:), delta_unscaled(:), atc(:), g(:), a(:,:), ata(:,:), L(:,:)
integer:: n, t, i, j

! Quantidade de estacoes e de parametros.
open(unit= 1, status= 'old', action= 'read', file= 'estacoes_e_parametros.dat')
read(1,*) n, t
close(1)

allocate (y(n), obs(n), previous_x(t), x(t), c(n), delta(t), delta_unscaled(t), a(n,t), atc(t), ata(t,t), g(n), L(t,t))

! Posicoes das estacoes.
open(unit= 2, status= 'old', action= 'read', file= 'estacoes.dat')
read(2,*) y(:)
close(2)

! Dados observados (sao dados sinteticos que foram gerados).
open(unit= 3, status= 'old', action= 'read', file= 'observado1.dat')
read(3,*) obs(:)
close(3)

! Valores iniciais dos parametros ("Chute").
open(unit= 4, status= 'old', action= 'read', file= 'kick_x1.dat')
read(4,*) previous_x(:)
close(4)

! Lambda
open(unit= 5, status= 'old', action= 'read', file= 'lambda1.dat')
read(5,*) lambda
close(5)

! Os dados sendo armazenados abaixo estao ai para observar o andamento do processo iterativo.
!---------------------------------------------------------------------------!
! Valores de x.
open(unit= 6, file='final_x.dat', status='old', action='write')

! Residuos.
open(unit= 7, file='c_residuo.dat', status='old', action='write')

! Delta.
open(unit= 8, file='delta.dat', status='old', action='write')

! Valores calculados.
open(unit= 9, file='calculado.dat', status='old', action='write')

! Matriz A.
open(unit= 10, file='a.dat', status='old', action='write')

! Matriz A.tc.
open(unit= 11, file='atc.dat', status='old', action='write')

! Merit function 's'.
open(unit= 12, file= 'merit_function.dat', status= 'old', action= 'write')

! Numero de iteracoes.
open(unit= 13, file= 'iteracoes.dat', status= 'old', action= 'write')

! Valores de lambda.
open(unit= 14, file= 'valores_lambda.dat', status= 'old', action= 'write')

! Valores de ata.
open(unit= 15, file= 'ata.dat', status= 'old', action= 'write')

!------------------------------------------------------------------------------!
previous_s= 0.
it= 0.
write(6,*) previous_x(:)
do  
    write(14,*) lambda
    it= 1 + it
    CALL gravity_and_c_values(n,t,previous_x,y,obs,g,c)
!------------------------------------------------------------------    
    ! Escreva os valores de g
    do i= 1, n
        write(9,*) g(i)
    end do
    
    ! Escreva os valores de c
    do i= 1, n
        write(7,*) c(i)
    end do    
!-----------------------------------------------------------------    
    ! Merit function. soma do quadrado dos resíduos.
    s= (sum(c**2))/(n-2)
    write(12,*) s   
    
    ! Calcula as matrizes ata e atc.
    CALL matrizes_ata_atc(n,t,lambda,previous_x,y,c,a,atc,ata)
!-------------------------------------------------------------------    
    ! Escreva a matriz a.
    do i= 1, n
        write(10,*) (a(i,j), j= 1, t)
    end do
    
    ! Escreva a matriz atc
    write(11,*) (atc(j), j= 1, t)
    
    ! Escreva a matriz ata
    do i= 1, t
        write(15,*) (ata(i,j), j= 1, t)
    end do
!--------------------------------------------------------------------
    ! Decompoe matriz ata.
    CALL sub1_decomp(t,ata,L)
    
    ! Soluciona o problema inverso para encontrar o novo valor de delta.
    CALL sub2_solve(t,ata,L,atc,delta)
    write(8,*) delta(:)  
    
    ! Eu não lembro por que eu criei uma matriz "delta_unscaled". Ela não dá o resultado esperado. Entao, vou deixar comentado até descobrir.
    ! Comentei esta iteração abaixo, mas ainda há memória alocada para esta matriz "delta_unscaled", lá no início do código.
    
    !do i= 1, t
    !    delta_unscaled(i)= delta(i)/(sqrt(ata(i,i)))
    !end do
    !write(8,*) delta_unscaled(:)   
    
    ! Parar a iteracao se:
    ! if ((abs(delta_unscaled(1)) <= +tol) .AND. (abs(delta_unscaled(2)) <= +tol)) exit
    if ((abs(delta(1)) <= +tol) .AND. (abs(delta(2)) <= +tol)) exit
    
    ! Se a condicao anterior nao for satisfeita:
    ! x(:)= previous_x(:) + delta(:)
    if (s < previous_s) then
        lambda= lambda/k
    else
        lambda= lambda
    end if
    ! x(:)= previous_x(:) + delta_unscaled(:)
    x(:)= previous_x(:) + delta(:)
    write(6,*) x(:)
    previous_x(:)= x(:)
    previous_s= s
    
end do

! Write out x final values.
write(6,*) '# valores estimados da profundidade(z) e do raio(a), em quilômetros (km), respectivamente:', (x(i), i= 1, t)
write(13,*) it
close(6)
close(7)
close(8)
close(9)
close(10)
close(11)
close(12)
close(13)
close(14)

deallocate (y,obs,previous_x,x,c,delta,delta_unscaled,a,ata,atc,l)

end program
