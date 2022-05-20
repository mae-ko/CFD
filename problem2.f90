dimension u(0:1000), flux(1000), x(1000)

!input
imax = 10
a = 1.
xmax = 1.
dx = xmax/float(imax)
cfl = 1

dt = cfl*dx/a

open(10, file = 'ini.d')

! numerical grod (why need
do i=1, imax 
  x(i)=dx*float(i)-0.5
enddo
!initial condition
do i = 1, imax-1
  xg = 0.5*(x(i) + x(i+1)) !why do this
 
  if (xg .le. 0.) then
     u(i) = 1.
  else
     u(i) = 0.
  endif
  !write(10, *) x(i),u(i)
 ! write(6,*) xg,u(i)
enddo

!nstep
do n = 1,5
!voundary condition
    u(0)=u(1)
    u(imax)=u(imax-1)

    do i = 1, imax
      flux(i) =0.5*a*(u(i-1)+u(i)) !central difference
      !flux(i) = a*(u(i-1)+u(i))/2 - abs(a)*(u(i)-u(i-1))/2  !upwond
    enddo
    
    do i = 1, imax-1
       u(i) = u(i) - dt/dx*(flux(i+1)-flux(i))
    enddo
    t = t + dt
enddo  

    x1 = a*t
    u1 = 1.
    x2 = a*t
    u2 = 0.
    write(6,*) x1,u1
    write(6,*) x2,u2
    write(6,*) 
   do i = 1, imax-1
    xg = 0.5*(x(i) + x(i+1))
    write(6,*) xg,u(i)
   enddo
!enddo
  
stop
end 
