! program: Hartree_Fock.f90
!
! This code performs self consistent HF calculation on a general lattice. The interaction
! term is decomposed in a rotationally symmetric fashion. 
!
! Author: Oinam Nganba Meetei
! Date:   02/15/12
!*********************************************************************************************************
module utilities
  use mtprng
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.d0)
  public::fermi, spiral
contains
  real(kind=double) function fermi(en,beta)
    implicit none
    real(kind=double),intent(IN)::en,beta
    fermi=1.d0/(exp(beta*en)+1.d0)
  end function fermi
 
  subroutine spiral(L,M,field)
    implicit none
    integer,intent(IN)::L,M
    complex(kind=double),dimension(:),intent(out)::field
    real(kind=double),dimension(:,:),allocatable::xyz
    integer::i,j,nsite
    real(kind=double)::phase
    open(unit=10,file='xyz.out',status='OLD')
    nsite=L*M
    allocate(xyz(nsite,3))
    do i=1,nsite
       read(10,*) (xyz(i,j),j=1,3)
    enddo
    field=0d0
    do i=1,2*nsite
       field(i)=0.5d0
    enddo
    do i=1,nsite
       phase=xyz(i,1)*4.d0*pi/3 + pi/5
       field(2*nsite+i)=0.5*exp((0d0,1d0)*phase)
    enddo
  end subroutine spiral
 
end module utilities



program hartree_fock
  use mtprng
  use hamiltonian
  use utilities
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.d0)
  real(kind=double)::dV,t,U,mu,beta,eta
  integer::seed
  character(40)::outfile

  integer:: nsite,i,j,count,INFO,ntarget,L,M
  real(kind=double)::err,tol,mu0,ntot,nup,ndn,dn,dmu,phase,energy,mixing
  
  real(kind=double),dimension(:),allocatable::eval,V,origfield
  complex(kind=double),dimension(:),allocatable::field,newfield
  real(kind=double),dimension(:,:),allocatable::ham
  complex(kind=double),dimension(:,:),allocatable::hmat
  real(kind=double),dimension(:,:),allocatable::xyz
  
  
  read(*,*) L,M
  read(*,*) dV,seed
  read(*,*) t, U
  read(*,*) mu,beta,eta
  read(*,*) outfile
  write(*,*) 't,U,dV',t,U,dV
  write(*,*) 'mu,beta',mu,beta
  write(*,*) 'eta',eta
  write(*,*) 'seed',seed
  write(*,*) 'filename  ',outfile


  ! Allocate variables
  nsite=L*M
  tol=1d-5
  mixing=0.3
  allocate(hmat(2*nsite,2*nsite),field(3*nsite),newfield(3*nsite),origfield(nsite), &
           eval(2*nsite),ham(nsite,nsite),V(nsite),xyz(nsite,3),stat=INFO)
  if(INFO.ne.0) then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif


  ! Generate non-interacting Hamiltonian for one spin sector
  call init_ham()
  call generate_disorder(dV,seed,V)
  mu0=mu
  mu=0d0
  call generate_hamiltonian(t,mu,V,ham)
  mu=mu0
  ntarget=nsite
  

  !TRIAL RUN FOR INITIALIZING MU
  !Initialize mean fields
  open(unit=1,file=outfile,status='OLD',iostat=INFO)
  if(info.eq.0) then
     do i=1,3*nsite
        read(1,*) field(i)
     enddo
     close(1)
     open(unit=1,file=outfile,status='REPLACE')
  else
     open(unit=1,file=outfile,status='NEW')
     open(unit=3,file='trial_'//outfile,status='OLD')
     call spiral(L,M,field)
     do i=1,nsite
        read(3,*) origfield(i)
     enddo
     field(1:nsite)=origfield
     field(nsite+1:2*nsite)=origfield
  endif

  
 
  ! Hartree-Fock loop
  count=0
  do
     !Generate full hamiltonian (contains both spins)
     hmat=0d0
     hmat(1:nsite,1:nsite)=ham
     hmat(nsite+1:2*nsite,nsite+1:2*nsite)=ham
     do i=1,2*nsite
        hmat(i,i)=hmat(i,i)+U*real(field(i))-mu        
     enddo
     do i=1,nsite
        hmat(i,nsite+i)=-U*field(2*nsite+i)
        hmat(nsite+i,i)=-U*conjg(field(2*nsite+i))
     enddo
     
     !Diagonalize full hamiltonian
     call diagonalize(hmat,eval)
     
     !Calculate new mean field parameters
     newfield=0d0
     do i=1,nsite
        do j=1,2*nsite
           !if(eval(j).gt.0) exit
           newfield(i)=newfield(i) + conjg(hmat(nsite+i,j))*hmat(nsite+i,j)*fermi(eval(j),beta)
           newfield(nsite+i)=newfield(nsite+i) + conjg(hmat(i,j))*hmat(i,j)*fermi(eval(j),beta)
        enddo
     enddo
     do i=1,nsite
        do j=1,2*nsite
           !if(eval(j).gt.0) exit
           newfield(2*nsite+i)=newfield(2*nsite+i) + conjg(hmat(nsite+i,j))*hmat(i,j)*fermi(eval(j),beta)
        enddo
     enddo
     
     
     !Calculate error in mean field parameters
     err=maxval(abs(field-newfield))
     write(*,*) 'err',err,count+1
     if(err.lt.tol) exit
     
     !update mean fields
     field=(1.d0-mixing)*newfield + mixing*field
     !field(1:nsite)=origfield
     !field(nsite+1:2*nsite)=origfield
     
     !Stop HF loop if it doesn't converge in 1000 steps
     count=count+1
     if (count.ge.20) goto 20
     
  enddo
  
  
20 ndn=real(sum(field(1:nsite)))
  nup=real(sum(field(nsite+1:2*nsite)))
  ntot=real(sum(field(1:2*nsite)))
  write(*,*) 'ntot,nup,ndn',ntot,nup,ndn
  
  do i=1,3*nsite
     write(1,*) field(i) 
  enddo

  !open(unit=3,file='spins_'//outfile)
  !do i=1,2*nsite
  !   write(3,*) real(field(i))
  !enddo
  !do i=1,nsite
  !   write(3,*) real(field(2*nsite+i))
  !enddo
  !do i=1,nsite
  !   write(3,*) imag(field(2*nsite+i))
  !enddo

  !Calculate MF energy
  energy=0d0
  do i=1,2*nsite
     energy=energy+eval(i)*fermi(eval(i),beta)
  enddo
  do i=1,nsite
     energy=energy-U*field(i)*field(nsite+i)+U*field(2*nsite+i)*conjg(field(2*nsite+i))
  enddo
  write(*,*) 'MF energy =', energy/nsite

  open(unit=20,file='eval'//outfile)
  do i=1,2*nsite
     write(20,*) i,eval(i)
  enddo

10 write(*,*) 'exiting HF'
  
end program Hartree_fock

