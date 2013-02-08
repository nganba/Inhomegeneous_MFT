program plot_neighbours
  
  implicit none
  integer :: i,j,k,n,n_neigh,job,junk
  real :: x,y,z
  integer, dimension(20)::neighbours

  open(unit=1,file='NN.out')
  open(unit=2,file='NNN.out')
  open(unit=3,file='xyz.out')
  open(unit=4,file='plot_neighbours.out')


  write(*,*) 'which point (enter an integer between 1 and L^3)'
  read(*,*) n
  
  do i=1,n
     read(3,*) x,y,z
  enddo
  write(4,*) x,y,z
  rewind(3)

  write(*,*) '1 for NN, 2 for NNN'
  read(*,*) job

  if(job.eq.1) then
     read(1,*) junk
     do i=1,n
        read(1,*) n_neigh
        read(1,*)(neighbours(j),j=1,n_neigh)
     enddo
     do i=1,n_neigh
        do j=1,neighbours(i)
           read(3,*) x,y,z
        enddo
        write(4,*) x,y,z
        rewind(3)
     enddo
  else
     read(2,*) junk
     do i=1,n
        read(2,*) n_neigh
        read(2,*)(neighbours(j),j=1,n_neigh)
     enddo
     do i=1,n_neigh
        do j=1,neighbours(i)
           read(3,*) x,y,z
        enddo
        write(4,*) x,y,z
        rewind(3)
     enddo
  endif

end program plot_neighbours
