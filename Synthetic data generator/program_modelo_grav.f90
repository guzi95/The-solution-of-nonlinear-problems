! Este codigo soluciona o seguinte problema direto: encontrar a componente gravitacional gerada 
! por uma esfera homogenea, em cada uma de 20 estacoes, espacadas entre si em 1 km. 
! A esfera tem raio 5 km, profundidade do centro igual a 7 km e densidade 0.25 g/cm3
! Autora: JULIANA DA SILVA FURTADO - BRASIL, PARÁ

program modelo_esfera_homogenea
implicit none
integer:: i, n
integer, parameter:: p= 2
real:: z, r    ! Profundidade e raio da esfera.
real, parameter:: pi= 3.1415, k= 6.672, d= 0.25    ! Valor de pi, constante gravitacional, contraste de densidade.
real, allocatable:: y(:), g(:)

! Le os dados matematicos do modelo (esfera homogenea).
open (file='parametros_a.dat', status='old', unit=10, action='read')
read(10,*) z, r, n    !Raio da esfera, profundidade do centro da esfera, numero de estacoes.
close(10)
allocate (y(n),g(n))

! Determina o espacamento de 1 km entre as estacoes.
open (file='estacoes.dat', unit= 11, action= 'write', status='old')
y(1)= -10
do i= 2, n
    y(i)= y(i-1)+1
enddo

do i= 1, 20
    write(11,*) y(i)
enddo

! A partir dos dados lidos anteriormente, obtem-se a anomalia gravitacional gerada pela esfera.
open (file= 'g_gerado.dat', unit= 12, status='new', action= 'write')
do i= 1, n
    g(i)= ((((4*pi)/3)*k*d*(r**3)*z)/((y(i)**2)+((z**2)**(3/2))))    ! O resultado final esta em mGal.
enddo
do i= 1, n
    write(12,100) y(i), g(i)
    100 format(' ', 10F12.6)
end do

deallocate (y,g)
endprogram
