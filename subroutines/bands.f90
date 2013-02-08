! Program:  spectral_function.f90
! 
!       This program calculates the single particle spectral function as a function of 
!       momentum and frequency for a given disorder realization. Disorder averaging is done
!       later
!
! Author: Oinam Nganba Meetei
! Date:   09/23/2011
!
!*************************************************************************************************
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



program spectral_function

  use common_func
  use hamiltonian
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  integer:: L,M,seed,nsite,i,j,nr,nk,n,nkpts
  real(kind=double)::U,dV,t,omega,eta,mu,beta,nup,qx,qy,x1,x2,y1,y2,tmp2,Akw
  complex(kind=double)::tmp
  real(kind=double),dimension(:,:),allocatable::xyz,ham,bz
  real(kind=double),dimension(:),allocatable::eval,V,pk
  complex(kind=double),dimension(:,:), allocatable::hmat
  complex(kind=double),dimension(:),allocatable::field
  character(50)::outfile


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
  

  open(unit=1,file='Band'//outfile,status='NEW')
  
  !Allocate variable
  nsite=L*M
  allocate(xyz(nsite,3),V(nsite),eval(2*nsite),field(3*nsite),ham(nsite,nsite),hmat(2*nsite,2*nsite))

  !Read mean-field parameters
  open(unit=2,file=outfile,status='OLD')
  read(2,*) mu
  do i=1,3*nsite
     read(2,*) field(i)
  enddo
  close(2)

  !Read co-ordinates
  open(unit=3,file='xyz.out',status='OLD')
  do i=1,nsite
     read(3,*)(xyz(i,j),j=1,3)
  enddo
  close(3)

  !Read k-points in 1st BZ
  open(unit=4,file='band_klist.out',status='OLD')
  i=0
  do
    read(4,*,end=10) tmp2
    i=i+1
  enddo
10 nkpts=i
  write(*,*) '# k points',nkpts
  allocate(bz(nkpts,3))
  rewind(4)
  do i=1,nkpts
     read(4,*) (bz(i,j),j=1,3)
  enddo
  close(4)
   

  !Generate non-interacting Hamiltonian
  call init_ham()
  call generate_disorder(dV,seed,V)
  call generate_hamiltonian(t,mu,V,ham)


  !Generate Mean-Field Hamiltonian
  hmat=0d0
  hmat(1:nsite,1:nsite)=ham
  hmat(nsite+1:2*nsite,nsite+1:2*nsite)=ham
  do i=1,2*nsite
     hmat(i,i)=hmat(i,i)+U*real(field(i))
  enddo
  do i=1,nsite
     hmat(i,nsite+i)=-U*field(2*nsite+i)
     hmat(nsite+i,i)=-U*conjg(field(2*nsite+i))
  enddo
 
 
  !Diagonalize Hamiltonian
  call diagonalize(hmat,eval)
    
  
  !CALCULATE SPECTRAL FUNCTION
  eta=0.1
  allocate(pk(2*nsite))
  do nk=1,nkpts
     qx=bz(nk,1); qy=bz(nk,2)
     do n=1,2*nsite
        tmp=0d0
        do i=1,nsite
           do j=1,nsite
              x1=xyz(i,1)
              y1=xyz(i,2)
              x2=xyz(j,1)
              y2=xyz(j,2)
              tmp = tmp + exp((0d0,-1.d0)*(qx*(x1-x2)+qy*(y1-y2)))*hmat(i,n)*conjg(hmat(j,n))
           enddo
        enddo
        pk(n) = real(tmp)
     enddo
     do omega=-10d0,10d0,0.1d0
        !write(1,*) 'omega',omega
        Akw=0d0
        do n=1,2*nsite
           Akw = Akw + pk(n)*lorentz(omega-eval(n),eta)/(2*nsite)
        enddo
        write(1,*) qx,omega,Akw
     enddo
  enddo

end program spectral_function
