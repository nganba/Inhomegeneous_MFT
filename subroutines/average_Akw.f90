
program average

  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer:: L,M,nsite,i,j,p,q,iter,INFO,nkpts
  real(kind=double)::omega,tmp
  character(40)::fname,infile
  character(5):: junk
  real(kind=double),dimension(:),allocatable::Akw,Akwtmp,kx,ky
 

  open(unit=1,file='names.dat',status='OLD')
  open(unit=11,file='avgMDCAkw.dat',status='NEW')


   i=0
   read(1,*) fname
   open(unit=2,file=fname,status='OLD')
   do
     read(2,*,end=1) tmp
     i=i+1
   enddo
1  nkpts=i
   rewind(1)
   close(2)
  
  allocate(Akw(nkpts),Akwtmp(nkpts),kx(nkpts),ky(nkpts),stat=INFO)
  if(INFO.ne.0)then
     write(*,*) 'Allocation of arrays failed'
     stop
  endif


   Akw=0d0
   iter=0
20 read(1,*,end=30) fname
   !write(*,*) fname 
   iter=iter+1
   open(unit=2,file=fname,status='OLD')
   do i=1,nkpts
      read(2,*,end=40) kx(i),ky(i),Akwtmp(i)
   enddo
   Akw=Akw+Akwtmp
   close(2)
   goto 20

30 Akw=Akw/iter
   write(*,*) '# copies averaged over',iter
   write(*,*) 'omega',omega
   !write(11,*) 'omega',omega
   do i=1,nkpts
      write(11,*) kx(i),ky(i),Akw(i)
   enddo
   write(11,*)

40 write(*,*) 'Averaging done'


end program average
