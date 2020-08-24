program Levenberg_Marquardt_Method
use levenberg_marq_mod
implicit none
real:: tol= 1.0E-5, k= 2, lambda, s, previous_s, it
real, allocatable:: y(:), obs(:), previous_x(:), x(:), c(:), delta(:), delta_unscaled(:), atc(:), g(:), a(:,:), ata(:,:), L(:,:)
integer:: n, t, i, j

! qnt de estacoes e de parametros
open(unit= 1, status= 'old', action= 'read', file= 'estacoes_e_parametros.dat')
read(1,*) n, t
close(1)

allocate (y(n),obs(n),previous_x(t),x(t),c(n),delta(t),delta_unscaled(t), a(n,t),atc(t),ata(t,t),g(n),L(t,t))

! posicoes das estacoes
open(unit= 2, status= 'old', action= 'read', file= 'estacoes.dat')
read(2,*) y(:)
close(2)

! dados observados (dados sinteticos)
open(unit= 3, status= 'old', action= 'read', file= 'observado1.dat')
read(3,*) obs(:)
close(3)

! Valores iniciais dos parametros ('chute').
open(unit= 4, status= 'old', action= 'read', file= 'initial_x1.dat')
read(4,*) previous_x(:)
close(4)

! lambda
open(unit= 5, status= 'old', action= 'read', file= 'lambdainicial.dat')
read(5,*) lambda
close(5)

! Os dados sendo armazenados abaixo estao ai para observar o andamento do processo iterativo.
!---------------------------------------------------------------------------!
! Valores de x.
open(unit= 6, file='final_x.dat', status='old', action='write')

! Resíduos.
open(unit= 7, file='c_resíduo.dat', status='old', action='write')

! Delta.
open(unit= 8, file='delta.dat', status='old', action='write')

! Valores calculados.
open(unit= 9, file='calculado.dat', status='old', action='write')

! Matriz a
open(unit= 10, file='a.dat', status='old', action='write')

! matriz atc
open(unit= 11, file='atc.dat', status='old', action='write')

! merit function 's'
open(unit= 12, file= 'merit_function.dat', status= 'old', action= 'write')

! numero iteracoes
open(unit= 13, file= 'iteracoes.dat', status= 'old', action= 'write')

! valores de lambda
open(unit= 14, file= 'valores_lambda.dat', status= 'old', action= 'write')

! valores de ata
open(unit= 15, file= 'ata.dat', status= 'old', action= 'write')

!------------------------------------------------------------------------------!
previous_s= 0.
it= 0.
write(6,*) previous_x(:)
do  
    write(14,*) lambda
    it= 1+it
    CALL gravity_and_c_values(n,t,previous_x,y,obs,g,c)
!------------------------------------------------------------------    
    ! escreva os valores de g
    do i= 1, n
        write(9,*) g(i)
    end do
    
    ! escreva os valores de c
    do i= 1, n
        write(7,*) c(i)
    end do    
!-----------------------------------------------------------------    
    ! Merit function. soma do quadrado dos resíduos.
    s= (sum(c**2))/(n-2)
    write(12,*) s   
    
    ! calcula as matrizes ata e atc.
    CALL matrizes_ata_atc(n,t,lambda,previous_x,y,c,a,atc,ata)
!-------------------------------------------------------------------    
    ! escreva a matriz a.
    do i= 1, n
        write(10,*) (a(i,j), j= 1, t)
    end do
    
    ! escreva a matriz atc
    write(11,*) (atc(j), j= 1, t)
    
    !escreva a matriz ata
    do i= 1, t
        write(15,*) (ata(i,j), j= 1, t)
    end do
!--------------------------------------------------------------------
    ! decompoe matriz ata
    CALL sub1_decomp(t,ata,L)
    
    ! soluciona o problema inverso para encontrar o novo valor de delta.
    CALL sub2_solve(t,ata,L,atc,delta)
    do i= 1, t
        delta_unscaled(i)= delta(i)/(sqrt(ata(i,i)))
    end do
    write(8,*) delta_unscaled(:)   
    
    ! Parar a iteracao se:
    if ((abs(delta_unscaled(1)) <= +tol) .AND. (abs(delta_unscaled(2)) <= +tol)) exit
    
    !se a condicao anterior nao for satisfeita:
    !x(:)= previous_x(:) + delta(:)
    
    if (s < previous_s) then
        lambda= lambda/k
    else
        lambda= lambda
    end if
    x(:)= previous_x(:) + delta_unscaled(:)
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
