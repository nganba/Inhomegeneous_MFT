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
  integer:: seed,nsite,i,j,i2,j2,p,q,l1,l2,iter,INFO,L,M,X,Y,col,indx,ix,iy
  real(kind=double)::U,dV,t,eta,mu,beta,mu0,norm,omega,avgLDOS,dx,dy,dr,r1,r2,LDOSofR
  character(40)::outfile
  real(kind=double),dimension(:,:),allocatable::ham,corr,corr2,xyz
  real(kind=double),dimension(:),allocatable::eval,V,LDOS
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

  !Generate and diagonalize hamiltonian
  call init_ham() 
  call generate_disorder(dV,seed,V)
  call generate_hamiltonian(t,mu,V,ham)
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
  write(*,*) 'Diagonalized'

  
  !Initialize local DOS arrays
  allocate(LDOS(2*nsite),stat=INFO)
  if(INFO.ne.0) then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif

  
  !Calculate local DOS
  LDOS=0d0
  omega=0d0  
  do i=1,2*nsite
     LDOS=LDOS + hmat(:,i)*conjg(hmat(:,i))*delta(eval(i)-omega,eta)
  enddo
  avgLDOS=sum(LDOS)/nsite
  do i=1,L
     do j=1,M
        write(10,*) i,j,LDOS((i-1)*M+j)+LDOS(nsite+(i-1)*M+j)-avgLDOS
     enddo
  enddo


  !Calculate correlation in LDOS
  allocate(corr(L,M))  !corr contains correlation of LDOS as a function 
                               !of relative distance
  

  corr=0d0
  do i=12,12!L
     do j=12,12!M
        do p=0,L-1
           do q=0,M-1
              l1=(i-1)*M+j
              i2=i+p
              j2=j+q
              if(i2.gt.L) i2=i2-L
              !if(i2.lt.1) i2=i2+L
              if(j2.gt.M) j2=j2-M
              !if(j2.lt.1) j2=j2+M
              l2=(i2-1)*M+j2
              ix=xyz(l2,1)-xyz(l1,1)
              iy=xyz(l2,2)-xyz(l1,2)              
              !corr(ix,iy)=corr(ix,iy) + (LDOS(l1)+LDOS(nsite+l1)-avgLDOS)*(LDOS(l2)+LDOS(nsite+l2)-avgLDOS)
              corr(p+1,q+1)=corr(p+1,q+1)+(LDOS(l1)+LDOS(nsite+l1)-avgLDOS)*(LDOS(l2)+LDOS(nsite+l2)-avgLDOS)
              write(50,*) ix,iy,(LDOS(l1)+LDOS(nsite+l1)-avgLDOS)*(LDOS(l2)+LDOS(nsite+l2)-avgLDOS)
           enddo
        enddo
     enddo 
  enddo
  !corr=corr/nsite
  do i=1,L
     do j=1,M
        !write(*,*) i,j,xyz((i-1)*M+j,1),xyz((i-1)*M+j,2)
        write(30,*) i,j,corr(i,j)
     enddo
  enddo

  !Perform angular averaging
  allocate(corr2(nsite,2))
  do i=1,L
     do j=1,M
        indx=(i-1)*M + j
        dx=i-1
        if(i-1 .gt. L/2) dx=L+1-i
        dy=j-1
        if(j-1 .gt. M/2) dy=M+1-j
        dr=sqrt(dx**2 + dy**2)
        corr2(indx,2)=corr(i,j)
        corr2(indx,1)=dr
     enddo
  enddo
  X=nsite;Y=2;col=1
  call sortby(corr2,X,Y,col)
  LDOSofR=0d0
  iter=0
  r1=corr2(1,1)
  do i=1,nsite
     r2=corr2(i,1)
     if(abs(r1-r2).lt.1d-4) then
        LDOSofR=LDOSofR+corr2(i,2)
        iter=iter+1
     else 
        write(40,*) r1,LDOSofR/iter
        LDOSofR=corr2(i,2)
        iter=1
        r1=r2
     endif
  enddo
  write(40,*) r1,LDOSofR/iter


end program DOSautocorrelation

