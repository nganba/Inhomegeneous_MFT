! generate_lattice.f90 
!
!      This program generates the lattice points of a crystal defined by the
! lattice vectors (to be read from input file). The size of lattice is fixed by 
! cubic cutoffs set (to be provided again in the input file) 
!
! Original: 02/09/2010  by Oinam Nganba Meetei
!
! Modifications:
!     02/09/12:: Generalized to any lattice

program generate_lattice

  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer :: i,j,k,Lx,Ly,Lz
  real(kind=double) :: x,y,z             
  real(kind=double),dimension(3,3):: basis   !a,b,c are the lattice parameters

  open(unit=2,file='xyz.out')

  read(*,*) Lx,Ly,Lz
  do i=1,3
     read(*,*)(basis(i,j),j=1,3)  !basis= co-ordinates for basis sites
  enddo
  write(*,*) 'Lx,Ly,Lz',Lx,Ly,Lz
  write(*,*) 'primitive vectors'
  do i=1,3
     write(*,*)(basis(i,j),j=1,3)
  enddo

  do i=0,Lx-1
     do j=0,Ly-1
        do k=0,Lz-1
           x = real(i)*basis(1,1) + real(j)*basis(2,1) + real(k)*basis(3,1)
           y = real(i)*basis(1,2) + real(j)*basis(2,2) + real(k)*basis(3,2)
           z = real(i)*basis(1,3) + real(j)*basis(2,3) + real(k)*basis(3,3)
           write(2,*) x,y,z
        enddo
     enddo
  enddo




  
  

!****************************************************************
  !generate for cubic lattice
  !do i=0,L-1
  !   do j=0,L-1
  !      do k=0,L-1
  !         do m=1,nbasis
  !            x=real(i)+basis(m,1)
  !            y=real(j)+basis(m,2)
  !            z=real(k)+basis(m,3)
  !            write(2,*) x,y,z     ! B  sublattice
  !            write(2,*) x+0.5,y,z ! B' sublattice
  !         enddo
  !       enddo
  !    enddo
  !enddo
!****************************************************************** 


end program generate_lattice
  

