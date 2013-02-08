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
  public::fermi, spiral, randomspin,neel
contains
  real(kind=double) function fermi(en,beta)
    implicit none
    real(kind=double),intent(IN)::en,beta
    fermi=1.d0/(exp(beta*en)+1.d0)
  end function fermi

  subroutine randomspin(L,M,field)
    implicit none
    integer,intent(IN)::L,M
    complex(kind=double),dimension(:),intent(out)::field
    integer::i,j,nsite
    real(kind=double)::phase,rand
    nsite=L*M
    field=0d0
    do i=1,2*nsite
       field(i)=0.5d0
    enddo
    do i=1,nsite
       phase=rand(0)
       field(2*nsite+i)=0.05*exp((0d0,1d0)*phase)
    enddo
  end subroutine randomspin

  subroutine neel(L,M,field)
    implicit none
    integer,intent(IN)::L,M
    complex(kind=double),dimension(:),intent(out)::field
    real(kind=double),dimension(:,:),allocatable::xyz
    integer::i,j,k,nsite
    real(kind=double)::phase
    open(unit=10,file='xyz.out',status='OLD')
    nsite=L*M
    allocate(xyz(nsite,3))
    do i=1,nsite
       read(10,*) (xyz(i,j),j=1,3)
    enddo
    field=0d0
    do i=1,nsite
       j=nint(xyz(i,1)); k=nint(xyz(i,2))
       field(i)=0.5d0 + 0.3d0*((-1d0)**(j+k))
       field(nsite+i)=1.d0-field(i)
    enddo    
    field(2*nsite+1:3*nsite)=0.001
  end subroutine neel 


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
       field(2*nsite+i)=0.05*exp((0d0,1d0)*phase)
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

  integer:: nsite,i,j,count,INFO,ntarget,L,M,MFloop
  real(kind=double)::err,tol,mu0,ntot,nup,ndn,dn,dmu,phase,energy,mixing,muvar,mu1,sx,sy,sz
  
  real(kind=double),dimension(:),allocatable::eval,V,origfield,mutable
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
  mixing=0.2
  MFloop=1500
  allocate(hmat(2*nsite,2*nsite),field(3*nsite),newfield(3*nsite),origfield(nsite), &
           eval(2*nsite),ham(nsite,nsite),V(nsite),xyz(nsite,3),mutable(MFloop),stat=INFO)
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
  ntarget=nsite !nint(nsite*1.01d0)
  

  !TRIAL RUN FOR INITIALIZING MU
  !Initialize mean fields
  open(unit=1,file=outfile,status='OLD',iostat=INFO)
  if(INFO.eq.0) then
     read(1,*) mu
     do i=1,3*nsite
        read(1,*) field(i)
     enddo
     close(1)
     open(unit=1,file=outfile,status='REPLACE')
  else
     !call spiral(L,M,field)
     !call randomspin(L,M,field)
     call neel(L,M,field)
     open(unit=1,file=outfile,status='NEW')
  endif
  
   
  mutable=0d0 
  ! Hartree-Fock loop
  count=0
  do
     !Generate full hamiltonian (contains both spins)
     hmat=0d0
     hmat(1:nsite,1:nsite)=ham
     hmat(nsite+1:2*nsite,nsite+1:2*nsite)=ham
     !Confine spins in xy plane
     !do i=1,nsite
     !   ntot=field(i)+field(i+nsite)
     !   field(i)=ntot/2.d0; field(i+nsite)=ntot/2.d0
     !enddo
     !Ensure total spin is zero
     !sx=0d0; sy=0d0; sz=0d0
     !do i=1,nsite
     !   sx=sx+real(field(2*nsite+i))
     !   sy=sy-imag(field(2*nsite+i))
     !   sz=sz+real(field(nsite+i)-field(i))/2.d0
     !enddo
     !write(*,*) sx,sy,sz
     !sx=sx/nsite; sy=sy/nsite; sz=sz/nsite
     !do i=1,nsite
     !   field(nsite+i)=field(nsite+i)-sz
     !   field(i)=field(i)+sz
     !   field(2*nsite+i)=field(2*nsite+i)-sx+(0d0,1.d0)*sy
     !enddo
     !sx=0d0; sy=0d0; sz=0d0
     !do i=1,nsite
     !   sx=sx+real(field(2*nsite+i))
     !   sy=sy-imag(field(2*nsite+i))
     !   sz=sz+real(field(nsite+i)-field(i))/2.d0
     !enddo
     !write(*,*) sx,sy,sz
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
     
     !update mu
     mutable(count+1)=mu
     if(int(count/20)*20.eq.count) mu = mu + (eval(ntarget)+eval(ntarget+1))/2.d0
     
     !Stop HF loop if it doesn't converge in 1000 steps
     count=count+1
     if (count.ge. MFloop) goto 20
     
  enddo
  
  
20 ndn=real(sum(field(1:nsite)))
  nup=real(sum(field(nsite+1:2*nsite)))
  ntot=real(sum(field(1:2*nsite)))
  write(*,*) 'ntot,nup,ndn',ntot,nup,ndn
  write(*,*) 'mu',mu
  !Calculate net spin
  sx=0d0; sy=0d0; sz=0d0
  do i=1,nsite
     sx=sx+real(field(2*nsite+i))
     sy=sy-imag(field(2*nsite+i))
     sz=sz+real(field(nsite+i)-field(i))/2.d0
  enddo
  write(*,*) "Net spin =",sqrt(sx**2 + sy**2 + sz**2)/nsite
  ! calculate variation in mu
  !mu=sum(mutable)/(count+1)
  !do i=1,count+1
  !   muvar=muvar + (mutable(i)-mu)**2
  !enddo
  !write(*,*) 'muvar',sqrt(muvar)/(count+1)
  !write(*,*) 'muvar',abs(mu-mu0)
  if(err.lt.tol) write(*,*) 'converged'


  !Write MF in output file
  write(1,*) mu
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

  !open(unit=20,file='eval'//outfile)
  do i=1,2*nsite
     write(20,*) i,eval(i)
  enddo


10 write(*,*) 'exiting HF'
  
end program Hartree_fock

