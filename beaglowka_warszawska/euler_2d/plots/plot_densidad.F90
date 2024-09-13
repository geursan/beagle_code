program escribe_plot
implicit none
character(len = 40) :: filein1,filepart1, coma
integer :: it, lt1,n_outputs
real :: tiempo, t_fin

n_outputs = 99
t_fin = 0.5
filein1  = 'out-000.dat'


open(unit = 10, file = 'grafica.p', status = 'unknown')
write(10,*) 'set xrange [0:1]'
write(10,*) 'set yrange [0:1]'
write(10,*) 'set cbrange [0.1:1.0]'
write(10,*) 'set title "density"'
write(10,*) 'set view map'
write(10,*) 'unset surface'
write(10,*) 'set pm3d at b'
write(10,*) 'set size square'
write(10,*) 'set terminal gif animate delay 5'
write(10,*) 'set output "densidad.gif"'
write(10,*) 'set style data lines'
write(10,*) 'set xlabel "x" font "Helvetica,15"'
write(10,*) 'set ylabel "y" font "Helvetica, 15"'


do it = 0, n_outputs
   lt1 = len(trim(filein1))
   call namechange(4,6,it,filepart1)
   filein1(lt1-6:lt1-4) = filepart1(1:3)


tiempo = real(it)*t_fin/real(n_outputs) 

write(10,*) 'splot ', '"',trim(filein1),'"',' u 1:2:3 title "time =',&
& tiempo,' s"'
end do


close(10)

end program escribe_plot


!=========================================================================
 subroutine namechange(lt1,lt2,it,fi)
!=========================================================================

   implicit none

   integer  , intent(in)  :: lt1, lt2, it
   character(len=10), intent(out) :: fi
   integer                   :: i, ns

   fi = ''

   ns = it
   do i =  lt2 - lt1 + 1, 1, -1
     fi(i:i) = char(ichar('0') + mod(ns, 10))
     ns      = ns / 10
   end do

 end subroutine namechange

!=========================================================================

