
! This program calculates the inverse problem by the Steepest-Descent Method. 
! From the data observed in the field (here are gravitational components) we want to
! estimate the depth (parameter 1) and the radius (parameter 2) of the sphere in subsurface through an iterative process.

program program_steepest_descent
use steepest_descent_method
implicit none
real, parameter:: tol= 1.0E-5
real:: u= 0., s= 0.
real, allocatable:: y(:), obs(:), initial_x(:), x(:), c(:), grad(:), delta(:), atc(:), g(:), a(:,:)
integer:: n, t, i, j   

! Number of stations, parameter and constant u.
open (unit=1, file='estacoes_e_parametros.dat', action='read', status='old')
 read(1,*) n, t
close(1)

! mi value.
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

! Valores iniciais dos parametros("chute").
open(unit= 4, action= 'read', file= 'initial_x1.dat', status= 'old')
read(4,*) initial_x(:)
close(4)

! valores finais dos parametros.
open(unit= 5, file='final_x.dat', status='old', action='write')

! Residuos.
open(unit= 6, file='c_resíduo.dat', status='old', action='write')

! gradiente.
open(unit= 7, file='grad.dat', status='old', action='write')

! Delta.
open(unit= 8, file='delta.dat', status='old', action='write')

! valores de gravidade calculados
open(unit= 9, file='calculado.dat', status='old', action='write')

! A matrix.
open(unit= 10, file='a.dat', status='old', action='write')

! At.C matrix.
open(unit= 11, file='atc.dat', status='old', action='write')

! Write out s(merit function).
open(unit= 12, file= 'merit_function.dat', status= 'old', action= 'write')

! Start estimate of the depth and radius parameters.

do
    ! Calculates the gravitational field (g) obtained from the parameter values, and the error (c) between these values and the observed data.
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
    
    ! Calculates the gradient (grad) and the step (delta) that will be given toward the minimum of the function.
    CALL atc_matrix(n,t,u,initial_x,y,c,grad,delta,a,atc)
    
    ! Write out A matrix.
    do i= 1, n
        write(10,*) (a(i,j), j= 1, t)
    end do
    
    ! Write out At.c matrix.
    write(11,*) (atc(j), j= 1, t)
    
    ! Write out gradient values.
    write(7,*) (grad(i), i= 1, t)
    
    ! Write out delta values.
    write(8,*) (delta(i), i= 1, t)
    
    ! x(p+1)
    x(:)= initial_x(:) + delta(:)   
    write(5,*) (x(i), i= 1, t)
    
    ! Condition to stop iteration:
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
