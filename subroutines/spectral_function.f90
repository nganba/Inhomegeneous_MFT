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
  real(kind=double),parameter::pi=acos(-1.d0)
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
  use k_points
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.0)
  integer:: L,M,seed,nsite,i,j,nr,nk,n,nkpts
  real(kind=double)::U,dV,t,omega,eta,mu,beta,nup,qx,qy,x1,x2,y1,y2,tmp2,Akw,dis,b,theta,phi,phasediv
  complex(kind=double)::tmp
  real(kind=double),dimension(:,:),allocatable::xyz,bz
  real(kind=double),dimension(:),allocatable::eval,V
  complex(kind=double),dimension(:,:), allocatable::hmat,ham
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
  

  open(unit=1,file='Akw'//outfile,status='NEW')
  
  !Allocate variable
  open(unit=2,file='NN.out',status='OLD')
  read(2,*) nsite
  close(2)
  allocate(xyz(nsite,3),V(nsite),eval(2*nsite),field(3*nsite),ham(nsite,nsite),hmat(2*nsite,2*nsite))

  !Read mean-field parameters
  if(abs(U).gt.1d-4) then
     open(unit=2,file=outfile,status='OLD')
     read(2,*) mu
     do i=1,3*nsite
        read(2,*) field(i)
     enddo
     close(2)
  else
     field=0d0
  endif

  !Read co-ordinates
  open(unit=3,file='xyz.out',status='OLD')
  do i=1,nsite
     read(3,*)(xyz(i,j),j=1,3)
  enddo
  close(3)
   
  !Initialize k-vectors
  allocate(bz(2*nsite,3))
  call init_kvec()
  
  !Initialize single particle Hamiltonian
  call init_ham()
  call generate_disorder(dV,seed,V)


  ! Loop over twisted phase
  phasediv=2
  do theta=0d0,2.d0*pi-pi/phasediv,2.d0*pi/phasediv
  do phi=0d0,2.d0*pi-pi/phasediv,2.d0*pi/phasediv
     write(*,*) 'theta,phi', theta,phi
     !Generate non-interacting Hamiltonian
     call generate_hamiltonian(t,mu,V,ham,theta,phi)
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
      
     !Generate k-points
     call generate_kpoints(theta,phi,bz,nkpts)

     !CALCULATE SPECTRAL FUNCTION
     !eta=1d-2
     omega=0d0
     do nk=1,nkpts
        Akw=0d0
        qx=bz(nk,1)
        qy=bz(nk,2)
        do n=1,2*nsite
           if(lorentz(omega-eval(n),eta).lt.0.001d0) cycle
           tmp=0d0
           do i=1,nsite
              do j=1,nsite
                 x1=xyz(i,1)
                 y1=xyz(i,2)
                 x2=xyz(j,1)
                 y2=xyz(j,2)
                 tmp = tmp + exp((0d0,-1.d0)*(qx*(x1-x2)+qy*(y1-y2)))*( hmat(i,n)*conjg(hmat(j,n)) + hmat(nsite+i,n)*conjg(hmat(nsite+j,n)) )
              enddo
           enddo
           Akw = Akw + real(tmp)*lorentz(omega-eval(n),eta)/nsite
        enddo
        write(1,*) qx,qy,Akw
     enddo 

  enddo
  enddo
      

end program spectral_function
