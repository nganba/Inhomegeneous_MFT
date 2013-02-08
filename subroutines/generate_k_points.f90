! Program: generate_k_points.f90
! 
! From the given primitive vectors in real space, this code generates the 
! k-points inside the first Brillouin zone.
!
! Author: Oinam Nganba Meetei
! Date:  02/10/12
!
!**************************************************************************

program k_points

  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.d0)
  integer :: i,j,k,Lx,Ly,Lz,n,iter
  real(kind=double) kx,ky,kz,a,b,dis,dis2
  real(kind=double),dimension(3)::a1,a2,a3,b1,b2,b3
  real(kind=double),dimension(20,3)::NN
  
  read(*,*) Lx,Ly,Lz
  read(*,*)(a1(j),j=1,3)
  read(*,*)(a2(j),j=1,3)
  read(*,*)(a3(j),j=1,3)
  write(*,*) 'Primitive vectors'
  write(*,*) (a1(i),i=1,3)
  write(*,*) (a2(i),i=1,3)
  write(*,*) (a3(i),i=1,3)
  
  a=sqrt(a1(1)**2 + a1(2)**2 + a1(3)**2)

  ! Calculate reciprocal lattice vectors
  b1(1)=a2(2)*a3(3)-a2(3)*a3(2)
  b1(2)=a2(3)*a3(1)-a2(1)*a3(3)
  b1(3)=a2(1)*a3(2)-a2(2)*a3(1)
  b1=((2*pi)/(a*dot_product(a1,b1)))*b1

  b2(1)=a3(2)*a1(3)-a3(3)*a1(2)
  b2(2)=a3(3)*a1(1)-a3(1)*a1(3)
  b2(3)=a3(1)*a1(2)-a3(2)*a1(1)
  b2=((2*pi)/(a*dot_product(a2,b2)))*b2

  b3(1)=a1(2)*a2(3)-a1(3)*a2(2)
  b3(2)=a1(3)*a2(1)-a1(1)*a2(3)
  b3(3)=a1(1)*a2(2)-a1(2)*a2(1)
  b3=((2*pi)/(a*dot_product(a3,b3)))*b3

  b=sqrt(b1(1)**2 + b1(2)**2 + b1(3)**2)

  write(*,*)
  write(*,*) 'Reciprocal lattice vectors'
  write(*,*) (b1(i),i=1,3)
  write(*,*) (b2(i),i=1,3)
  write(*,*) (b3(i),i=1,3)


 
  !Find nearest neighbors in reciprocal lattice
  !Assumes that only on length 'a' defines NN sites uniquely
  NN=0d0
  n=0
  do i=-1,1
     do j=-1,1
        do k=-1,1
           kx=i*b1(1) + j*b2(1) + k*b3(1)
           ky=i*b1(2) + j*b2(2) + k*b3(2)
           kz=i*b1(3) + j*b2(3) + k*b3(3)
           dis=sqrt(kx**2 + ky**2 + kz**2)
           if(abs(dis-b).lt.1d-2) then
              n=n+1
              NN(n,1)=kx
              NN(n,2)=ky
              NN(n,3)=kz
           endif
        enddo
     enddo
  enddo
  do i=1,n
     write(20,*)(NN(i,j),j=1,3)
  enddo

  ! Sample kpoints
  open(unit=2,file='brillouin_zone.out')
  do i=-Lx,Lx
     do j=-Ly,Ly
        inner:do k=-Lz,Lz
           kx=i*b1(1)/Lx + j*b2(1)/Ly + k*b3(1)/Lz
           ky=i*b1(2)/Lx + j*b2(2)/Ly + k*b3(2)/Lz
           kz=i*b1(3)/Lx + j*b2(3)/Ly + k*b3(3)/Lz
           dis=sqrt(kx**2 + ky**2 + kz**2)
           !Check if k belongs to 1st BZ
           if(abs(kz).gt.1d-4) cycle inner ! use this for 2D lattice
           do iter=1,n
              dis2=sqrt((kx-NN(iter,1))**2 + (ky-NN(iter,2))**2 + (kz-NN(iter,3))**2)
              if(dis2.lt.dis) cycle inner
           enddo
           write(2,*)kx,ky,kz
        enddo inner
     enddo
  enddo
  
end program k_points
  
  
