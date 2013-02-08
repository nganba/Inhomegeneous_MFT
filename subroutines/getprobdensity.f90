! program: getprobden.f90
!
! This code produces the density profile for states close to w=0
!
! Author: Oinam Nganba Meetei
! Date:   07/24/12
!*********************************************************************************************************
program getprobden
  use hamiltonian
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  integer:: seed,nsite,i,j,count,INFO,ntarget,L,M,izero
  real(kind=double)::U,dV,t,eta,mu,beta,dos,omega,mu0
  character(40)::outfile
  real(kind=double),dimension(:,:),allocatable::ham,xyz
  real(kind=double),dimension(:),allocatable::eval,V
  complex(kind=double),dimension(:,:),allocatable::hmat
  complex(kind=double),dimension(:),allocatable::field
  
  
  read(*,*) L,M
  read(*,*) dV,seed
  read(*,*) t, U
  read(*,*) mu,beta,eta
  read(*,*) outfile
  write(*,*) L,M
  write(*,*) dV,seed
  write(*,*) t, U
  write(*,*) mu,beta, eta
  write(*,*) outfile
  
  nsite=L*M

  ! Allocate variables
  allocate(eval(2*nsite),hmat(2*nsite,2*nsite),ham(nsite,nsite),V(nsite),field(3*nsite),xyz(nsite,3),stat=INFO)
  if(INFO.ne.0) then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif

  !Read mean-field parameters
  if(abs(U).gt.1d-4) then
     write(*,*) 'Reading mean field'
     open(unit=1,file=outfile,status='OLD')
     read(1,*) mu
     do i=1,3*nsite
        read(1,*) field(i)
     enddo
     close(1)
  else
     field=0d0; mu=0d0
  endif

  !Read co-ordinates
  open(unit=2,file='xyz.out',status='OLD')
  do i=1,nsite
     read(2,*)(xyz(i,j),j=1,3)
  enddo
  close(2)

  ! Generate non-interacting Hamiltonian 
  call init_ham()
  call generate_disorder(dV,seed,V)
  call generate_hamiltonian(t,mu,V,ham)
  
  
  ! Generate and diagonalize HF Hamiltonian 
  hmat=0d0
  hmat(1:nsite,1:nsite)=ham
  hmat(nsite+1:2*nsite,nsite+1:2*nsite)=ham
  do i=1,2*nsite
     hmat(i,i)=hmat(i,i) + U*real(field(i))
  enddo
  do i=1,nsite
     hmat(i,nsite+i)=-U*field(2*nsite+i)
     hmat(nsite+i,i)=-U*conjg(field(2*nsite+i))
  enddo
  call diagonalize(hmat,eval)
  do i=1,2*nsite
     write(10,*) eval(i)
  enddo
  write(*,*) 'mu, (E(n)+E(n+1))/2',mu,(eval(nsite)+eval(nsite+1))/2.d0
  
  !Use this only for getting the non-interacting prob density
  if(abs(U).lt.1d-4) then
     mu=(eval(nsite)+eval(nsite+1))/2.d0
     eval=eval-mu
  endif
 

  izero=2*nsite
  do i=1,2*nsite
     if(abs(eval(i)).lt.abs(eval(izero))) izero=i
  enddo
  do j=1,nsite
     write(20,*) xyz(j,1),xyz(j,2), abs(hmat(j,izero))**2 + abs(hmat(j+nsite,izero))**2
     write(30,*) xyz(j,1),xyz(j,2), V(j)
     write(40,*) xyz(j,1),xyz(j,2), V(j) + U*real(field(j))
  enddo
  write(*,*) 'eigenvalue',eval(izero)
  
end program getprobden

