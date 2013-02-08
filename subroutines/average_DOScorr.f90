
program average

  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer:: i,j,L,iter,INFO
  character(40)::fname,infile
  character(5):: junk
  real(kind=double),dimension(:),allocatable::omega,dostmp,dos,r


  open(unit=1,file='names.dat',status='OLD')
  open(unit=11,file='avgDOScorr.dat',status='NEW')

   !Determine size of omega array
   read(1,*) fname
   write(*,*) fname
   open(unit=2,file=fname,status='OLD')
   i=0
   do
     read(2,*,end=10) junk,junk,junk
     i=i+1
   enddo
10 L=i
   close(2)
   rewind(1)
   allocate(dos(L),dostmp(L),omega(L),r(L),stat=INFO)
   if(INFO.ne.0)then
      write(*,*) 'Allocation of arrays failed'
      stop
   endif


   dos=0d0
   iter=0
20 read(1,*,end=30) fname
   !write(*,*) fname 
   iter=iter+1
   open(unit=2,file=fname,status='OLD')
   do i=1,L
     read(2,*) omega(i),r(i),dostmp(i)
   enddo
   dos=dos+dostmp
   close(2)
   goto 20

30 dos=dos/iter
   do i=1,L
      write(11,*) omega(i),r(i),dos(i)
   enddo


   write(*,*) 'Averaging done'


end program average
