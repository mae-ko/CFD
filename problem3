dimension u(0:1000),flux(1000),x(1000),c2(1000),r2(1000),r1(1000),u_delta(1000)
real g(100,100)
!input
imax = 100
a = 1.
xmax = 1.
dx = xmax/float(imax)
cfl = 1.

dt = cfl*dx/a

open(10, file = 'ini.d')

! numerical grod (why need
do i=1, imax 
  x(i) = dx*float(i)-0.5
enddo
!initial condition
do i = 1, imax-1
  xg = 0.5*(x(i) + x(i+1)) !why do this
 
  if (xg .le. 0.) then
     u(i) = 1.
  else
     u(i) = 0.
  endif
enddo

!nstep
do n = 1,15
!boundary condition
    u(0)=u(1)
    u(imax)=u(imax-1)

    do i = 1, imax
      flux(i) =0.5*a*(u(i-1)+u(i)) !central difference
      !flux(i) = a*(u(i-1)+u(i))/2 - abs(a)*(u(i)-u(i-1))/2  !upwond
    enddo
    
    do i = 1, imax-1
       r1(i) = -(flux(i+1)-flux(i))*dt/dx
    enddo
call thomas(u,g,c2,r2,r1,u_delta)
    do i = 1, imax-1
       u(i) = u(i) + u_delta(i)
    enddo
    t = t+dt
enddo  

    x1 = a*t
    u1 = 1.
    x2 = a*t
    u2 = 0.
    write(6,*) x1,u1
    write(6,*) x2,u2
    write(6,*) 
   do i = 1,imax-1
      xg = 0.5*(x(i) + x(i+1))
    write(6,*) xg,u(i)
   enddo
  
stop
end 

subroutine thomas (u,g,c2,r2,r1,u_delta)
dimension u(0:1000),c2(1000),r2(1000),u_delta(1000),r1(1000)
real g(100,100)
imax = 100
a = 1.
xmax = 1.
dx = xmax/float(imax)
cfl = 1.

dt = cfl*dx/a
!matrix
do k = 1,imax-1
   do j = 1,imax-1
   if (j==k) then
      g(j,j) = 1.
   else if (j==k+1) then
      g(j,k) = -0.5*dt/dx
   else if (j==k-1) then
      g(j,k) = 0.5*dt/dx
   else
      g(j,k) = 0.
   endif
   enddo
enddo

!boundary condition
g(1,1) = 1-0.5*dt/dx
g(imax-1,imax-1) = 1+0.5*dt/dx

!main
  c2(1) = g(1,2)/g(1,1)
  r2(1) = r1(1)/g(1,1)

  do i = 2, imax-2
     c2(i) = g(i,i+1)/(g(i,i)-c2(i-1)*g(i,i-1))
  enddo

  do i = 2, imax-1
     r2(i) = (r1(i)-r2(i-1)*g(i,i-1))/(g(i,i)-c2(i-1)*g(i,i-1))
  enddo
  
      u_delta(imax-1) = r2(imax-1)

  do i = imax-2, 1, -1
      call exit
      u_delta(i) = r2(i)-c2(i)*u_delta(i+1)
  enddo
end subroutine
