! Program:  brillouin_zone.f90
! 
!       This program calculates the k-points inside the first brillouin zone. 
!       
!
! Author: Oinam Nganba Meetei
! Date:   08/29/2012
!
!*************************************************************************************************



program brillouin_zone 

  use k_points
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  integer:: nsite,i,j,nkpts,L,M
  real(kind=double),dimension(:,:),allocatable::bz,basis
  real(kind=double)::theta,phi
  

  !Allocate arrays
  open(unit=2,file='triangular_lat_param.inp',status='OLD')
  read(2,*) L,M
  nsite=L*M
  close(2)
  allocate(bz(2*nsite,3),basis(3,3))

  !Real space primitive basis vectors
  basis(1,1)=1.5d0; basis(1,2)=-0.866025; basis(1,3)=0d0
  basis(2,1)=1.5d0; basis(2,2)=0.866025;  basis(2,3)=0d0
  basis(3,1)=0d0;   basis(3,2)=0d0;       basis(3,3)=1.d0

  !Initialize k-vectors
  call init_kvec(basis)
  
  !Generate k-points in 1st Brillouin zone
  theta=0d0; phi=0d0
  call generate_kpoints(theta,phi,bz,nkpts)
  open(unit=2,file='brillouin_zone.out')
  do i=1,nkpts
     write(2,*)(bz(i,j),j=1,3)
  enddo

end program brillouin_zone
