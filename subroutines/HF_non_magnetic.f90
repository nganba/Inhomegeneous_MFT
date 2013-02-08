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
  public::fermi
contains
  real(kind=double) function fermi(en,beta)
    implicit none
    real(kind=double),intent(IN)::en,beta
    fermi=1.d0/(exp(beta*en)+1.d0)
  end function fermi
 
 end module utilities



program hartree_fock
  use mtprng
  use hamiltonian
  use utilities
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.d0)
  real(kind=double)::dV,t,U,mu,beta,eta,mixing
  integer::seed
  character(40)::outfile

  integer:: nsite,i,j,count,INFO,ntarget,L,M
  real(kind=double)::err,tol,mu0,ntot,dn,dmu,phase,energy
  
  real(kind=double),dimension(:),allocatable::eval,V
  real(kind=double),dimension(:),allocatable::field,newfield
  real(kind=double),dimension(:,:),allocatable::ham
  real(kind=double),dimension(:,:),allocatable::hmat
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
  mixing=0.5
  allocate(hmat(nsite,nsite),field(nsite),newfield(nsite), &
           eval(nsite),ham(nsite,nsite),V(nsite),xyz(nsite,3),stat=INFO)
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
  ntarget=nsite/2
  

  !TRIAL RUN FOR INITIALIZING MU
  !Initialize mean fields
  open(unit=1,file='trail_'//outfile,status='OLD',iostat=INFO)
  if(INFO.eq.0) then
     do i=1,nsite
        read(1,*) field(i)
     enddo
     close(1)
     open(unit=1,file='trial_'//outfile,status='REPLACE')
  else
     open(unit=1,file='trial_'//outfile,status='NEW')
     if(abs(dV).gt.1d-4) then
        field=0.5d0-0.5*V/(dV)
     else
        field=0.5d0
     endif
  endif
  !field=ntarget*field/(sum(field))
 
  ! Hartree-Fock loop
  count=0
  do
     hmat=ham
     do i=1,nsite
        hmat(i,i)=hmat(i,i)+U*field(i)-mu-0.1        
     enddo
     
     !Diagonalize full hamiltonian
     call diagonalize(hmat,eval)
     
     !Calculate new mean field parameters
     newfield=0d0
     do i=1,nsite
        do j=1,nsite
           !if(eval(j).gt.0) exit
           newfield(i)=newfield(i) + hmat(i,j)*hmat(i,j)*fermi(eval(j),beta)
        enddo
     enddo
     
     
     !Calculate error in mean field parameters
     err=maxval(abs(field-newfield))
     write(*,*) 'err',err,count+1
     if(err.lt.tol) exit
     
     !update mean fields
     field=(1.d0-mixing)*newfield + mixing*field


     !calculate new mu
     !if(count.gt.100) mu=(eval(ntarget+1)+eval(ntarget))/2.d0
     
     !Stop HF loop if it doesn't converge in 2000 steps
     count=count+1
     if (count.ge.200) goto 20
     
  enddo
  
  !Update chemical potential
  ntot=real(sum(field(1:nsite)))
  dn=ntot-ntarget
  if(abs(dn).lt.1d-4) goto 20
  !mu0=mu
  !mu=mu+1d-3*mtprng_rand_real3(state)
  !dmu=mu-mu0
  write(*,*) 'ntot,ntarget',ntot,ntarget
  
  
  !PROPER RUN WITH LOOP FOR FIXING MU
  !do 
     ! Initialize mean fields
     !   "field" is a vector with all the site dependent fields. Elements 1 to nsite are 
     !   Hartree fields for up spin. (nsite+1) to 2*nsite are Hartree fields for dn spin.
     !   (2*nsite+1) to 3*nsite are <S-> fields which is the same as <S+> fields
  !   field=0.5d0
 
     ! Hartree-Fock loop
  !   count=0
  !   do
  !      hmat=ham
  !      do i=1,nsite
  !         hmat(i,i)=hmat(i,i)+U*field(i)-mu        
  !      enddo
  !      
  !      !Diagonalize full hamiltonian
  !      call diagonalize(hmat,eval)
  !      
  !      !Calculate new mean field parameters
  !      newfield=0d0
  !      do i=1,nsite
  !         do j=1,nsite
  !            !if(eval(j).gt.0) exit
  !            newfield(i)=newfield(i) + hmat(i,j)*hmat(i,j)*fermi(eval(j),beta)
  !         enddo
  !      enddo
  !              
  !      
  !      !Calculate error in mean field parameters
  !      err=maxval(abs(field-newfield))
  !      write(*,*) 'err',err,count+1
  !      if(err.lt.tol) exit
  !      
  !      !update mean fields
  !      field=(1.d0-mixing)*newfield + mixing*field
  !            
  !      !Stop HF loop if it doesn't converge in 2000 steps
  !      count=count+1
  !      if (count.ge.2000) goto 10
  !
  !   enddo
  !   
  !   !Update chemical potential
  !   ntot=real(sum(field(1:nsite)))
  !   if(abs(ntot-ntarget).lt.1d-4) then
  !      exit
  !   else
  !      if(abs(ntot-ntarget-dn).lt.1d-4) then
  !         if(ntot.lt.ntarget) then
  !            ntot=0
  !         else
  !            ntot=nsite
  !         endif
  !      endif
  !      mu=mu-(ntot-ntarget)*dmu/(ntot-ntarget-dn)
  !      dmu=mu-mu0
  !      !if(abs(dmu).gt.10) mu = mu0 + 0.01*(ntarget-ntot)*mtprng_rand_real3(state)
  !      dmu=mu-mu0
  !      mu0=mu
  !      dn=ntot-ntarget
  !   endif
  !   write(*,*)'ntot,ntot-ntarget',ntot,dn
  !   write(*,*)'mu,dmu',mu,dmu
  !  
  !   
  !enddo

20 ntot=real(sum(field(1:nsite)))
  write(*,*) 'ntot,ntarget',ntot,ntarget
  !write(*,*) 'mu',mu
  
  do i=1,nsite
     write(1,*) field(i)
  enddo
  !open(unit=2,file='V'//outfile)
  !do i=1,nsite
  !   write(2,*) V(i)
  !enddo  


10 write(*,*) 'exiting HF'
  
end program Hartree_fock

