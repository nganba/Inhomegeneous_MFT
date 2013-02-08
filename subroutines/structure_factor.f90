!  This code calculates the magnetic structure factor from the HF
!  solution. 

program structure_factor
  use k_points
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer::i,j,k,L,M,nsite,nkpts,seed
  real(kind=double)::theta,phi,qx,qy,x1,x2,y1,y2
  real(kind=double)::dV,t,U,mu,beta,eta
  complex(kind=double)::SF
  real(kind=double),dimension(:,:),allocatable::xyz,bz
  real(kind=double),dimension(:),allocatable::sx,sy,sz
  complex(kind=double),dimension(:),allocatable::field
  character(50)::outfile
  
  read(*,*) dV,seed
  read(*,*) t, U
  read(*,*) mu,beta,eta
  read(*,*) outfile

  open(unit=1,file=outfile,status='OLD')
  open(unit=2,file='triangular_lat_param.inp',status='OLD')
  read(2,*) L,M
  nsite=L*M
  close(2)

  allocate(xyz(nsite,3),bz(2*nsite,3),sx(nsite),sy(nsite),sz(nsite),field(3*nsite))

  !Read mean field
  read(1,*) mu
  do i=1,3*nsite
     read(1,*) field(i)
  enddo
  do i=1,nsite
     sz(i)=real(field(nsite+i)-field(i))/2.d0
     sx(i)=real(field(2*nsite+i))
     sy(i)=-imag(field(2*nsite+i))
  enddo
  open(unit=10,file='spins.out')
  do i=1,nsite
     write(10,*) sx(i),sy(i),sz(i)
  enddo
  
  !Read real-space coordinates
  open(unit=2,file='xyz.out',status='OLD')
  do i=1,nsite
     read(2,*)(xyz(i,j),j=1,3)
  enddo
  close(2)

  
  !Generate k-points
  call init_kvec()
  theta=0d0; phi=0d0
  call generate_kpoints(theta,phi,bz,nkpts)

  !Calculate structure factor
  open(unit=3,file='SF'//outfile,status='NEW')
  do i=1,nkpts
     SF=0d0
     qx=bz(i,1); qy=bz(i,2)
     do j=1,nsite
        do k=1,nsite
           x1=xyz(j,1); y1=xyz(j,2)
           x2=xyz(k,1); y2=xyz(k,2)
           SF = SF + ( sx(j)*sx(k) + sy(j)*sy(k) + sz(j)*sz(k) )*exp( (0d0,-1d0)*( qx*(x1-x2) + qy*(y1-y2) ) )
        enddo
     enddo
     if(imag(SF).gt.1d-6) then
        write(*,*) 'SF is not real'
        stop
     endif
     write(3,*) qx,qy,real(SF)/(nsite*nsite)
  enddo


end program structure_factor
