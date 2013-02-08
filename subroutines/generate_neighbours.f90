! generate_neighbours.f90
!
!     This program reads the co-ordinates of lattice sites and generates a 
!     nearest neighbor table. Periodic boundary condition is used.
!
! Original: 02/09/2010 by Oinam Nganba Meetei
!
! Modified: 02/09/2012:: generalized to any lattice defined by lattice vector

program neighbours
  
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer :: i,j,m1,m2,m3,nmax,n_NN,Nx,Ny,Nz
  real(kind=double):: r_NN,dis
  integer, dimension(:,:),allocatable::NN
  real,dimension(:),allocatable :: x,y,z
  real(kind=double),dimension(3,3):: basis 
  real(kind=double),dimension(3)::r

  open(unit=2,file='xyz.out')
  open(unit=3,file='NN.out')
  open(unit=4,file='NNN.out')
 
  read(*,*) Nx,Ny,Nz
  do i=1,3
     read(*,*)(basis(i,j),j=1,3)  !basis= co-ordinates for basis sites
  enddo
  r_NN=sqrt(basis(1,1)**2 + basis(1,2)**2 + basis(1,3)**2)
  write(*,*) 'NN distance',r_NN
  nmax=0
  do 
     read(2,*,END=1) dis,dis,dis
     nmax=nmax+1
  enddo
1 rewind(2)
  allocate(x(nmax))
  allocate(y(nmax))
  allocate(z(nmax))
  allocate(NN(nmax,50))
  

  do i=1,nmax
     read(2,*,end=10) x(i),y(i),z(i)
  enddo

10 write(*,*) 'total number of lattice points',nmax 

  !do i=1,n
  !   write(*,*) x(i),y(i),z(i)
  !enddo

  NN=0
  do i=1,nmax
     n_NN=0

     do j=1,nmax        
        if (i.eq.j) cycle

        !Check for NN
PbNN:   do m1=-1,1
           do m2=-1,1
              do m3=-1,1
                 r(1) = x(i)-x(j)-m1*Nx*basis(1,1)-m2*Ny*basis(2,1)-m3*Nz*basis(3,1)
                 r(2) = y(i)-y(j)-m1*Nx*basis(1,2)-m2*Ny*basis(2,2)-m3*Nz*basis(3,2)
                 r(3) = z(i)-z(j)-m1*Nx*basis(1,3)-m2*Ny*basis(2,3)-m3*Nz*basis(3,3)
                 dis = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
                 if(abs(dis-r_NN).le.1d-2) then
                    n_NN = n_NN + 1
                    NN(i,1) = n_NN
                    NN(i,n_NN+1)=j
                    exit PbNN
                 endif
              enddo
           enddo
        enddo PbNN

     enddo
  !write(*,*) i,NN(i,1),NNN(i,1)
  enddo

  write(*,*) 'genertated'
  write(3,120) nmax
  write(4,120) nmax
  do i=1,nmax
     write(3,120) NN(i,1)
     write(3,110)(NN(i,j),j=2,NN(i,1)+1)
  enddo
110 format(20(I6,1x))
120 format(I6)
  write(*,*) 'written'

end program neighbours


