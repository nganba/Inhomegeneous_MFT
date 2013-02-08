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
module utilities
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  public:: delta
contains
  function delta(x,eta)
     implicit none
     real(kind=double),intent(in)::x,eta
     real(kind=double)::delta
     delta=(1.d0/pi)*eta/(x*x+eta*eta)
  end function delta
end module utilities


program DOSautocorrelation
  use hamiltonian
  use utilities
  use sort_util
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  integer:: seed,nsite,i,j,i2,j2,iter,INFO,L,M,col,phasediv
  real(kind=double)::U,dV,t,eta,mu,beta,mu0,norm,omega,avgLDOS,dx,dy,dr,r1,r2,LDOSofR
  real(kind=double)::theta,phi
  character(40)::outfile
  character(1)::opt
  real(kind=double),dimension(:,:),allocatable::corr,xyz
  real(kind=double),dimension(:),allocatable::eval,V,LDOS
  complex(kind=double),dimension(:,:),allocatable::hmat,ham
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
  allocate(eval(2*nsite),hmat(2*nsite,2*nsite),ham(nsite,nsite),V(nsite),field(3*nsite),xyz(nsite,3),stat=INFO)
  if(INFO.ne.0) then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif

  !Read positions
  open(unit=1,file='xyz.out',status='OLD')
  do i=1,nsite
     read(1,*)(xyz(i,j),j=1,3)
  enddo
  close(1)

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


  !Initialize local DOS arrays
  allocate(LDOS(2*nsite),corr(nsite*nsite,2),stat=INFO)
  if(INFO.ne.0) then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif



  open(unit=1,file='DOScorr'//outfile,status='NEW')
  !Initialize  hamiltonian
  call init_ham()
  call generate_disorder(dV,seed,V)

  !Calculate local DOS
  do omega=-1,1,0.1
  LDOS=0d0
  !omega=0d0
  ! Loop over twisted phase
  phasediv=1
  do theta=0d0,2.d0*pi-pi/phasediv,2.d0*pi/phasediv
  do phi=-0d0,2.d0*pi-pi/phasediv,2.d0*pi/phasediv
     call generate_hamiltonian(t,mu,V,ham,theta,phi)
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
     !write(*,*) 'Diagonalized'

     !Calculate local DOS
     do i=1,2*nsite
        LDOS=LDOS + hmat(:,i)*conjg(hmat(:,i))*delta(eval(i)-omega,eta)
     enddo
  enddo
  enddo
  avgLDOS=sum(LDOS)/nsite
  !do i=1,L
  !   do j=1,M
  !      write(10,*) i,j,LDOS((i-1)*M+j)+LDOS(nsite+(i-1)*M+j)-avgLDOS
  !   enddo
  !enddo

  !Calculate autocorrelation
  corr=0d0
  iter=0
  do i=1,nsite
     do j=1,nsite
        iter=iter+1
        dx=xyz(j,1)-xyz(i,1)
        if(abs(dx).gt.L-abs(dx)) dx=L-abs(dx)
        dy=xyz(j,2)-xyz(i,2)
        if(abs(dy).gt.M-abs(dy)) dy=M-abs(dy)
        dr=sqrt(dx**2 + dy**2)
        corr(iter,1)=dr
        corr(iter,2)=(LDOS(i)+LDOS(nsite+i)-avgLDOS)*(LDOS(j)+LDOS(nsite+j)-avgLDOS)
     enddo 
  enddo
  !do i=1,nsite*nsite
  !   write(20,*) corr(i,1),corr(i,2)
  !enddo 
  i2=nsite*nsite;j2=2;col=1
  opt='a'
  call sortby(corr,i2,j2,col,opt)
  !do i=1,nsite*nsite
  !   write(30,*) corr(i,1),corr(i,2)
  !enddo
  LDOSofR=0d0
  iter=0
  r1=corr(1,1)
  do i=1,nsite*nsite
     r2=corr(i,1)
     if(abs(r1-r2).lt.1d-4) then
        LDOSofR=LDOSofR+corr(i,2)
        iter=iter+1
     else 
        write(1,*) omega,r1,LDOSofR/iter
        LDOSofR=corr(i,2)
        iter=1
        r1=r2
     endif
  enddo
  write(1,*) omega,r1,LDOSofR/iter
  enddo


end program DOSautocorrelation

