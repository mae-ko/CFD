dimension u(0:1000), dif(1000)

!input
imax = 50
a = 1.
xmax = 1.
dx = xmax/float(imax)
cfl = 0.5

dt = cfl*dx/a
 
open(10, file = 'initial.d')

!initial condition 
do i = 1, imax-1
    x = dx*float(i) - 0.5
   if (x .le. -0.2) then
       u(i)=0
   else if (x .gt. -0.2 .and. x .le. 0) then
       u(i)=1.+ 5.*x
   else if (x .gt. 0. .and. x .le. 0.2) then
       u(i)=1.-5*x
   else
       u(i)=0.
   endif
   write(10,*) x,u(i)
   !write(6,*) x,u(i)
enddo
!call exit
do n = 1, 10
!boundary condition
u(0)=u(1)
u(imax-1)=u(imax)

    do i = 1, imax-1
       dif(i)=0.5*(u(i+1)-u(i-1))
    enddo
    do i = 1, imax-1
       u(i)=u(i)-dt/dx*dif(i)
    enddo
    t = t + dt
enddo

    x1 = a*t - 0.2
    u1 = 0
    x2 = a*t
    u2 = 1
    x3 = a*t + 0.2
    u3 = 0
    write(6,*) x1, u1
    write(6,*) x2, u2
    write(6,*) x3, u3
    write(6,*)


 do i=1, imax-1
    x = dx *float(i)-0.5
 write(6,*) x, u(i)
 enddo
stop
end
