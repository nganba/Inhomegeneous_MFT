! program: single_particle.f90
!
! This code calculates the single particle eigenstates and eigenvalues. It uses
! the subroutines in hamiltonian.f90 to generate single particle Hamiltonian on
! a general lattice with oniste disorder. 
!
! Author: Oinam Nganba Meetei
! Date:   02/10/12
!*********************************************************************************************************
module common_func
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  
  public:: lorentz

contains

  real(kind=double) function lorentz(omega,eta)
    implicit none
    real(kind=double),parameter::pi=acos(-1.0)
    real(kind=double),intent(IN)::omega,eta
    lorentz = eta/(omega**2 + eta**2)
    lorentz = lorentz/pi
  end function lorentz

end module common_func




program single_particle
  use hamiltonian
  use common_func
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  integer:: seed,nsite,i,j,count,INFO,ntarget,L,M,n
  real(kind=double)::U,dV,t,eta,mu,beta,dos,omega,qx,qy,x1,x2,y1,y2,Akw
  complex(kind=double)::tmp
  character(40)::outfile
  real(kind=double),dimension(:,:),allocatable::hmat,xyz
  real(kind=double),dimension(:),allocatable::eval,V,pk
  
  
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
  allocate(eval(nsite),hmat(nsite,nsite),xyz(nsite,3),V(nsite),stat=INFO)
  if(INFO.ne.0) then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif
  

  !Read coordinates
  open(unit=2,file='xyz.out',status='OLD')
  do i=1,nsite
     read(2,*) (xyz(i,j),j=1,3)
  enddo
  close(2)


  ! Generate and diagonalize Hamiltonian 
  call init_ham()
  call generate_disorder(dV,seed,V)
  mu=0d0
  call generate_hamiltonian(t,mu,V,hmat)
  call diagonalize(hmat,eval)
  
  mu=(eval(nsite/2) + eval((nsite/2)+1))/2.d0
  eval=eval-mu
  open(unit=2,file="eval.dat")
  do i=1,nsite
     write(2,*) i,eval(i)
  enddo
  close(2)

  ! Calculate DOS
  !open(unit=3,file='sDOS'//outfile)
  !eta=eta
  !do omega=-10.d0,10.d0,0.02d0
  !   dos=0d0
  !   do i=1,nsite
  !      dos=dos+ (1.d0/pi)*eta/((omega-eval(i))**2 + eta**2)
  !   enddo
  !   write(3,*) omega,dos/nsite
  !enddo

  !Calculate spectral function
  open(unit=1,file='sAkw'//outfile)
  !Calculate eigenstate projection
  allocate(pk(nsite))
  qx=0d0; qy=2.8214669932
  do n=1,nsite
     tmp=0d0
     do i=1,nsite
        do j=1,nsite
           x1=xyz(i,1)
           y1=xyz(i,2)
           x2=xyz(j,1)
           y2=xyz(j,2)
           tmp = tmp + exp((0d0,-1.d0)*(qx*(x1-x2)+qy*(y1-y2)))*hmat(i,n)*hmat(j,n)
        enddo
     enddo
     pk(n) = real(tmp)
  enddo
  eta=1d-2
  do omega=-10d0,10d0,0.01d0
     !write(1,*) 'omega',omega
     Akw=0d0
     do n=1,nsite
        Akw = Akw + pk(n)*lorentz(omega-eval(n),eta)/(2*nsite)
     enddo
     write(1,*) omega,Akw
  enddo  


  !Print the prob density of state at Fermi energy
  !ntarget=nsite/2  
  !do j=1,nsite
  !   write(20,*) xyz(j,1),xyz(j,2), abs(hmat(j,ntarget))**2
  !   write(30,*) xyz(j,1),xyz(j,2), V(j) 
  !enddo  

  !open(unit=4,file='disorder.out')
  !do i=1,nsite
  !   write(4,*) i,V(i)
  !enddo
  ! Calculate dispersion
  
  
end program Single_particle

