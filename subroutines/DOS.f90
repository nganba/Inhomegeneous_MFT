! program: single_particle.f90
!
! This code calculates the DOS of HF Hamiltonian 
!
! Author: Oinam Nganba Meetei
! Date:   03/09/12
! 
! Modified: 08/15/12
!         Calculates DOS with twisted boundary conditions to get more eigenvalues.
!*********************************************************************************************************

program density_of_states
  use hamiltonian
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  integer:: seed,nsite,i,j,count,INFO,ntarget,L,M,nw,phasediv!,iter
  real(kind=double)::U,dV,t,eta,mu,beta,mu0,theta,phi,wmin,wmax,dw,norm
  character(40)::outfile
  complex(kind=double),dimension(:,:),allocatable::ham
  real(kind=double),dimension(:),allocatable::eval,V,omega,DOS!,evaltot
  complex(kind=double),dimension(:,:),allocatable::hmat
  complex(kind=double),dimension(:),allocatable::field
  
  
  read(*,*) L,M
  read(*,*) dV,seed
  read(*,*) t, U
  read(*,*) mu,beta,eta
  read(*,*) outfile
  write(*,*) dV,seed
  write(*,*) t, U
  write(*,*) mu,beta, eta
  write(*,*) outfile
  
  nsite=L*M

  ! Allocate variables
  allocate(eval(2*nsite),hmat(2*nsite,2*nsite),ham(nsite,nsite),V(nsite),field(3*nsite),stat=INFO)
  if(INFO.ne.0) then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif

  !Read mean-field parameters
  if(abs(U).gt.1d-4) then
     open(unit=1,file=outfile,status='OLD')
     read(1,*) mu
     do i=1,3*nsite
        read(1,*) field(i)
     enddo
     close(1)
  else
     field=0d0
  endif

  !Initialize DOS arrays
  dw=0.005d0
  wmax=5.d0
  wmin=-5.d0
  nw=nint((wmax-wmin)/dw) + 1
  allocate(omega(nw),DOS(nw))
  do i=1,nw
     omega(i)= wmin + (i-1)*dw
  enddo
  DOS=0d0

  !Initialize single particle hamiltonian
  call init_ham() 
  call generate_disorder(dV,seed,V)

  ! Loop over twisted phase
  phasediv=4
  !allocate(evaltot(2*nsite*((2*phasediv)**2)))
  !iter=0
  do theta=0d0,2.d0*pi-pi/phasediv,2.d0*pi/phasediv
  do phi=-0d0,2.d0*pi-pi/phasediv,2.d0*pi/phasediv
     write(*,*) 'theta,phi', theta,phi
     ! Generate non-interacting Hamiltonian 
     mu0=mu
     mu=0d0
     call generate_hamiltonian(t,mu,V,ham,theta,phi)
     mu=mu0
    
     ! Generate and diagonalize HF Hamiltonian 
     hmat=0d0
     hmat(1:nsite,1:nsite)=ham
     hmat(nsite+1:2*nsite,nsite+1:2*nsite)=ham
     do i=1,2*nsite
        hmat(i,i)=hmat(i,i) + U*real(field(i)) - mu
     enddo
     do i=1,nsite
        hmat(i,nsite+i)=-U*field(2*nsite+i)
        hmat(nsite+i,i)=-U*conjg(field(2*nsite+i))
     enddo
     call diagonalize(hmat,eval)
     !eta=2d-2
     do i=1,nw
        do j=1,2*nsite
           DOS(i)=DOS(i)+ (1.d0/pi)*eta/((omega(i)-eval(j))**2 + eta**2)
        enddo
     enddo
     !do i=1,2*nsite
     !   evaltot( iter*2*nsite + i)=eval(i)
     !enddo
     !iter=iter+1
  enddo
  enddo

  do i=1,2*nsite
     write(10,*) i,eval(i)
  enddo
  ! Calculate DOS
  !eta=eta
  !do i=1,nw
  !   do j=1,2*phasediv*2*nsite
  !      DOS(i)=DOS(i)+ (1.d0/pi)*eta/((omega(i)-evaltot(j))**2 + eta**2)
  !   enddo
  !enddo
 
  !Normalize DOS to 1
  norm=0d0
  do i=1,nw
     norm = norm + DOS(i)*dw 
  enddo
  DOS=DOS/norm

  !write DOS in file
  open(unit=3,file='DOS'//outfile)
  do i=1,nw
     write(3,*) omega(i),DOS(i)
  enddo

end program density_of_states

