module my_params
    double precision           :: c         = 0.9d0  !courant number
    integer, parameter         :: n_outputs = 100 !n outputs
    double precision           :: tfin      = 0.5d0  !tiempo final
    integer, parameter         :: nx        = 100    !n celdas
    double precision, parameter:: Lx        = 1d0    !tamano del dominio
    double precision :: deltaX = Lx/nx
    double precision :: gam = 1.4d0  !gamma del gas
!
end module
!
program euler_1d
use my_params
    implicit none
    double precision :: deltaT,t
    integer :: i, count_output
    double precision :: up(5,nx)
    double precision :: u(5,0:nx+1), uL(5,0:nx+1), uR(5,0:nx+1) ! vectores u(x)
    double precision :: F(3,0:nx+1) !vectores
    double precision :: x
    !
    !definir condiciones iniciales, ie, llenar el vector u inicial
    !recuerda !ee = P / (rho * (gam - 1d0) )
              !E  = rho*(0.5d0*v**2 + ee)
    do i=1,nx
        x=real(i)*deltaX   ! la posicion en x
        !un escalon
        if(x.le.0.3d0) then
           u(4,i) = 0.75d0    !velocidad
           u(5,i) = 1d0       !presion
           u(1,i) = 1d0       !densidad
           u(2,i) = u(1,i)*u(4,i)
           u(3,i) = u(1,i)*(0.5d0*u(4,i)**2 + u(5,i) / ( u(1,i) * (gam - 1d0) ) ) 
        else 
           u(4,i) = 0d0         !velocidad
           u(5,i) = 0.1d0       !presion
           u(1,i) = 0.125d0     !densidad
           u(2,i) = u(1,i)*u(4,i)
           u(3,i) = u(1,i)*(0.5d0*u(4,i)**2 + u(5,i) / ( u(1,i) * (gam - 1d0) ) )  
        end if
    
        if(u(1,i).lt.0d0.or.u(5,i).lt.0d0)then
        print*, 'presiones o densidades negativas en las condiciones iniciales'
        print*, 'presion:',u(5,i), 'densidad:', u(1,i)
        stop
        end if
    end do
!
!     
    t = 0d0          !tiempo
    i = 0            !contador
    count_output = 0 !contador para imprimir archivos
    ! evolucion de los vectores
    do while(t.le.tfin) !loop para evolucion temporal
      
      u(:,0) = u(:,1)
      !
      u(:,nx+1) = u(:,nx)  
      !
      ! buscar sL sR

      sL = min(vxL-csL,vxR-csR)
      sR = max(vxL-csL,vxR-csR)

      !Definir UL y UR
      call prims(sL,sR,u,uL,uR)  
      ! Flujos iniciales
      call fluxes(u,F)
      !
      !paso de tiempo
      call timestep(u,i,deltaT)
      !
      !integrador      
      call Lax_integrator(t,deltaT,u,F,count_output,up)
      ! 
      u(:,1:nx)  = up(:,:)  !actualiza la "u", todo el vector
      !
      t = t + deltaT
      !
    end do
    !
end program euler_1d
!
!===========================================================================
subroutine prims(u,uL,uR)
use my_params 
implicit none
double precision, intent(in)  :: u(5,0:nx+1), sL, sR, time 
double precision, intent(out) :: uL(5,0:nx+1), uR(5,0:nx+1) ! vectores u(x)
double precision :: x
integer :: i
uL=0d0
uR=0d0
do i = 0, nx +1
x=real(i)*deltaX 
 if(x.le.(sR*time)) then
 uL(:,i)=u(:,i)
 else if(x.ge.(sL*time))
 uR(:,i)=u(:,i)
end do
end do

end subroutine
!===========================================================================

!===========================================================================
subroutine fluxes(u,F)
use my_params
implicit none
integer :: j
double precision, intent(in)  :: u(5,0:nx+1)
double precision, intent(out) :: F(3,0:nx+1)



  do j=0,nx+1
  !
    F(1,j) = u(1,j)*u(4,j)
	F(2,j) = u(1,j)*u(4,j)**2 + u(5,j)
    F(3,j) = u(4,j)*( u(3,j) + u(5,j) )
  !
  end do
end subroutine
!===========================================================================

!===========================================================================
subroutine timestep(u,i,deltaT)
use my_params
implicit none
integer, intent(inout) :: i
double precision, intent(in) :: u(5,0:nx+1)
double precision :: cs(0:nx+1)
double precision :: deltaT,t
integer :: j
!velocidad del sonido en cada celda
  !
  i = i + 1   !contador
  !
  do j = 0, nx+1
  !
    cs(j) = sqrt(gam*u(5,j) / u(1,j))
    !   
    if(u(1,j).lt.0d0.or.u(5,j).lt.0d0)then
    print*, 'presiones o densidades negativas'
    print*, u(5,j), u(1,j), 'i=', i
    stop
    end if
  !
  end do
  !
  !
  ! Paso de tiempo
  deltaT = c*deltaX/maxval(abs(u(4,:))+cs(:))
  !
  if(i.le.5)then
  deltaT=0.2d0*deltaT
  end if
  !
end subroutine
!===========================================================================

!===========================================================================
subroutine Lax_integrator(t,deltaT,u,F,count_output,up)
!esta subroutina es el integrador e imprime el archivo
use my_params
implicit none
 double precision, intent(in) :: t, deltaT, F(3,0:nx+1), u(5,0:nx+1)
 integer, intent(inout) :: count_output 
 double precision, intent(out):: up(5,nx)
 double precision :: x
 integer :: j, k
 character (len=20) output
!
  if(t.ge.(count_output*tfin/real(n_outputs)))then
  !
  write(output,'(a,i3.3,a)') 'out-',count_output,'.dat'
  count_output = count_output + 1
  open(unit=10,file=output,status='unknown')
  end if
  !
  !
  !imprimir un vector "u" en cada paso de tiempo
  !
   do j=1,nx !loop para llenar espacios X
     x=real(j)*deltaX   ! la posicion en x
     !
     do k = 1, 3
     up(k,j) = 0.5d0*( u(k,j-1) + u(k,j+1) - ( F(k,j+1)-F(k,j-1) )*deltaT/deltaX )  !LAX
     end do
     !
     up(4,j) = up(2,j) / up(1,j) !velocidad (momento/rho)
     up(5,j) = ( up(3,j) - 0.5d0*up(1,j)*up(4,j)**2 ) * (gam - 1d0) !presion
     ! 
     !
     !posision, densidad, energia, velocidad, presion
     write(10,*) x, up(1,j), up(3,j), up(4,j), up(5,j)
         ! 
   end do
   close(10)
end subroutine
!===========================================================================
