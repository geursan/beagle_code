module my_params
    double precision           :: c         = 0.4d0 ! 0.9d0  !courant number
    integer, parameter         :: n_outputs = 100    !n outputs
    double precision           :: tfin      = 0.5d0  !tiempo final
    integer, parameter         :: nx        = 200    !n celdas
    integer, parameter         :: ny        = 200    !n celdas
    double precision, parameter:: Lx        = 1d0    !tamano del dominio
    double precision, parameter:: Ly        = 1d0    !tamano del dominio
    double precision :: deltaX = Lx/nx
    double precision :: deltaY = Ly/ny
    double precision :: gam = 1.4d0  !gamma del gas
!
end module
!
program euler_2d
use my_params
    implicit none
    double precision :: deltaT,t
    integer :: i,j,count_output
    double precision :: u(7,0:nx+1,0:ny+1) !vectores u -->
    !             !u(label de cantidad conservada, valores en x, valores en y )
    double precision :: up(7,nx,ny) !para actualizar conservadas a partir de flujos.
    !             !up(label de cantidad conservada, valores en x, valores en y )
    double precision :: F(4,0:nx+1,0:ny+1) !Flujos F
    double precision :: G(4,0:nx+1,0:ny+1) !Flujos G
    double precision :: x, y, r, vx, vy, vr
    !
    !definir condiciones iniciales, ie, llenar el vector u inicial
    !recuerda !ee = P / (rho * (gam - 1d0) )
              !E  = rho*(0.5d0*v**2 + ee)
    vx=0.75d0
    vy=0.75d0
    vr = sqrt(vx**2+vy**2)
    ! condiciones iniciales en las celdas
    !---------------------------------------------------------------------------
    do i=1,nx
        x=real(i)*deltaX   ! la posicion en x
      do j=1,ny
        y=real(j)*deltaY   ! la posicion en y
        !un escalon
        r = sqrt(x**2+y**2)
        !if(r.le.0.3d0.and.r.gt.0d0.and.(atan(abs(x/y)).le.0.26d0)) then
        if(r.le.0.3d0.and.r.gt.0d0.and.(atan(abs(x/y)).ge.(0.78d0-0.26d0)).and.(atan(abs(x/y)).le.(0.78d0+0.26d0))) then
        !if(r.le.0.3d0.and.r.gt.0d0) then
           u(5,i,j) = vr*x/r    !velocidad en x
           u(6,i,j) = vr*y/r    !velocidad en y
           u(7,i,j) = 1d0 !1d0      !presion
           u(1,i,j) = 1d0 !1d0      !densidad
           u(2,i,j) = u(1,i,j)*u(5,i,j)
           u(3,i,j) = u(1,i,j)*u(6,i,j)
           u(4,i,j) = u(1,i,j)*(0.5d0*(u(5,i,j)**2+u(6,i,j)**2) + u(7,i,j) / ( u(1,i,j) * (gam - 1d0) ) )
        else
           u(5,i,j) = 0d0    !velocidad en x
           u(6,i,j) = 0d0    !velocidad en y
           u(7,i,j) = 0.1d0    !presion
           u(1,i,j) = 0.125d0  !densidad
           u(2,i,j) = u(1,i,j)*u(5,i,j)
           u(3,i,j) = u(1,i,j)*u(6,i,j)
           u(4,i,j) = u(1,i,j)*(0.5d0*(u(5,i,j)**2+u(6,i,j)**2) + u(7,i,j) / ( u(1,i,j) * (gam - 1d0) ) )
        end if

        if(u(1,i,j).lt.0d0.or.u(7,i,j).lt.0d0)then
        print*, 'presiones o densidades negativas en las condiciones iniciales'
        print*, 'presion:',u(7,i,j), 'densidad:', u(1,i,j)
        stop
        end if
      end do
    end do
    !---------------------------------------------------------------------------
!
!
    t = 0d0          !tiempo
    i = 0            !contador
    count_output = 0 !contador para imprimir archivos
    ! evolucion de los vectores
    do while(t.le.tfin) !loop para evolucion temporal
      !
     
      u(:,1:nx,0)    = u(:,1:nx,1)  !actualiza fila debajo
      u(:,1:nx,ny+1) = u(:,1:nx,ny) !actualiza fila arriba
      u(:,0,1:ny)    = u(:,1,1:ny)  !columna izquierda
      u(:,nx+1,1:ny) = u(:,nx,1:ny) !columna derecha
      !
      u(:,0,0)       = u(:,0,1)     !actualiza una esquina de la caja
      u(:,nx+1,0)    = u(:,nx,0)    !actualiza una esquina de la caja
      u(:,0,ny+1)    = u(:,1,ny+1)  !actualiza una esquina de la caja
      u(:,nx+1,ny+1) = u(:,nx,ny+1) !actualiza una esquina de la caja
      !
      ! Flujos iniciales
      call fluxes(u,F,G)
      !
      !paso de tiempo
      call timestep(u,i,deltaT)
      !
      !integrador
      call Lax_integrator(t,deltaT,u,F,G,count_output,up)
      !
      u(:,1:nx,1:ny)  = up(:,:,:)  !actualiza la "u", todo el vector
      !
      t = t + deltaT
      !
    end do
    !
end program euler_2d
!
!===========================================================================
subroutine fluxes(u,F,G)
use my_params
implicit none
integer :: j,k
double precision, intent(in)  :: u(7,0:nx+1,0:ny+1)
double precision, intent(out) :: F(4,0:nx+1,0:ny+1), G(4,0:nx+1,0:ny+1)
  do j=0,nx+1
    do k=0,ny+1
    !
    F(1,j,k) = u(1,j,k)*u(5,j,k)
    F(2,j,k) = u(1,j,k)*u(5,j,k)**2 + u(7,j,k)
	F(3,j,k) = u(1,j,k)*u(5,j,k)*u(6,j,k)
    F(4,j,k) = u(5,j,k)*( u(4,j,k) + u(7,j,k) )
    !
    G(1,j,k) = u(1,j,k)*u(6,j,k)
    G(2,j,k) = u(1,j,k)*u(5,j,k)*u(6,j,k)
    G(3,j,k) = u(1,j,k)*u(6,j,k)**2 + u(7,j,k)
    G(4,j,k) = u(6,j,k)*( u(4,j,k) + u(7,j,k) )
    !
    end do
  end do
end subroutine
!===========================================================================

!===========================================================================
subroutine timestep(u,i,deltaT)
use my_params
implicit none
integer, intent(inout) :: i
double precision, intent(in) :: u(7,0:nx+1,0:ny+1)
double precision :: cs(0:nx+1,0:ny+1)
double precision :: deltaT,t,valx,valy
integer :: j, k
!velocidad del sonido en cada celda
  !
  i = i + 1   !contador
  !
  do j = 0, nx+1
    do k = 0, ny+1
    !
    cs(j,k) = sqrt(gam*u(7,j,k) / u(1,j,k))
    !
    if(u(1,j,k).lt.0d0.or.u(7,j,k).lt.0d0)then
    print*, 'presiones o densidades negativas'
    print*,u(7,j,k), u(1,j,k), 'i=', i
    stop
    end if
    !
    end do
  end do
  !
  !
  valx = c*deltaX/maxval(abs(u(5,:,:))+cs(:,:))
  valy = c*deltaY/maxval(abs(u(6,:,:))+cs(:,:))
  ! Paso de tiempo
  deltaT = min(valx,valy)
  !
  if(i.le.5)then
  deltaT=0.2d0*deltaT
  end if
  !
end subroutine
!===========================================================================

!===========================================================================
subroutine Lax_integrator(t,deltaT,u,F,G,count_output,up)
!esta subroutina es el integrador e imprime el archivo
use my_params
implicit none
 double precision, intent(in) :: t, deltaT, u(7,0:nx+1,0:ny+1)
 double precision, intent(in) :: F(4,0:nx+1,0:ny+1), G(4,0:nx+1,0:ny+1)
 integer, intent(inout) :: count_output
 double precision, intent(out):: up(7,nx,ny)
 double precision :: x, y
 integer :: j, k, l
 character (len=20) output
!
  if(t.ge.(count_output*tfin/real(n_outputs)))then
  !
  print*, 'output=',count_output
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
     do k=1,ny
      y=real(k)*deltaY   ! la posicion en x
      !
      do l = 1, 4
        !LAX
      up(l,j,k) = 0.25d0*( u(l,j+1,k) + u(l,j,k+1) + u(l,j-1,k) + u(l,j,k-1) ) &
                & -0.5d0*( F(l,j+1,k) - F(l,j-1,k) )*deltaT/deltaX  &
                & -0.5d0*( G(l,j,k+1) - G(l,j,k-1) )*deltaT/deltaY
      end do
      !
      up(5,j,k) = up(2,j,k) / up(1,j,k) !velocidad x (momento en x /rho)
      up(6,j,k) = up(3,j,k) / up(1,j,k) !velocidad y (momento en y /rho)
      up(7,j,k) = ( up(4,j,k) - 0.5d0*up(1,j,k)&
                  &*(up(5,j,k)**2 + up(6,j,k)**2) ) * (gam - 1d0) !presion
      !
      !
      !posision x,posision y,densidad,energia,velocidad x,velocidad y,presion
      !write(10,*) x, y, up(1,j,k), up(4,j,k), up(5,j,k), up(6,j,k), up(7,j,k)
      write(10,'(6es13.5)') x, y, up(1,j,k), up(4,j,k), up(5,j,k), up(6,j,k), up(7,j,k)
            !
     end do
     write(10,*)
   end do
   close(10)
end subroutine
!===========================================================================
