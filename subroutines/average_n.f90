module utilities
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.d0)
  public::fermi
contains
  real(kind=double) function fermi(en,beta)
    implicit none
    real(kind=double),intent(IN)::en,beta
    fermi=1.d0/(exp(beta*en)+1.d0)
  end function fermi
end module utilities


program average_n
  use utilities
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer:: i,j,L,M,iter,INFO,seed,nstate
  real(kind=double)::dV,t,U,mu,beta,eta,tmpn,avgn
  character(40)::fname,infile,outfile
  real(kind=double),dimension(:),allocatable::eval


  open(unit=1,file='names.dat',status='OLD')

  !Read data
  read(*,*) L,M
  read(*,*) dV,seed
  read(*,*) t, U
  read(*,*) mu,beta,eta
  read(*,*) outfile
 
  nstate=2*L*M
  allocate(eval(nstate))

  avgn=0d0
  iter=0
20 read(1,*,end=30) fname
  !write(*,*) fname 
  iter=iter+1
  open(unit=2,file=fname,status='OLD')
  do i=1,nstate
     read(2,*) j,eval(i)
  enddo
  tmpn=0d0
  do i=1,nstate
     tmpn=tmpn+fermi(eval(i),beta)
  enddo
  write(*,*) tmpn
  avgn=avgn+tmpn
  close(2)
  goto 20

30 avgn=avgn/iter
   write(*,*) 'avg n =', avgn


   write(*,*) 'Averaging done'


 end program average_n
