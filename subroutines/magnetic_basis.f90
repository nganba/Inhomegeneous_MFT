! triangular_magnetic_basis.f90 
!
!      This program identifies the non-equivalent sites on a triangular lattice
!      with 120 degree spin order
!
! Author: Oinam Nganba Meetei
! Date  : 08/29/12

module magnetic_basis

  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer,save::Lx,Ly,Lz,nsite
  real(kind=double),dimension(:,:),allocatable,save::xyz

  public::identify_magnetic_basis

contains

  subroutine identify_magnetic_basis(mag_basis)
    implicit none
    integer,dimension(:),intent(out)::mag_basis
    integer :: i,j,id             
    real(kind=double),dimension(3,3):: basis,d   !a,b,c are the lattice parameters
    real(kind=double),dimension(3)::r1,r2
 
    !Read lattice parameters and initialize nearest neighbor vectors
    open(unit=2,file='triangular_lat_param.inp',status='OLD')
    read(2,*) Lx,Ly,Lz
    do i=1,3
       read(2,*)(basis(i,j),j=1,3)  !basis= primitive basis vectors
    enddo
    close(2)
    d(1,:)=basis(1,:)
    d(2,:)=basis(2,:)
    d(3,:)=basis(2,:)-basis(1,:)
  
    !Read co-ordinates of real space lattice
    nsite=Lx*Ly*Lz
    allocate(xyz(nsite,3))
    open(unit=2,file='xyz.out',status='OLD')
    do i=1,nsite
       read(2,*)(xyz(i,j),j=1,3)
    enddo
    close(2)
    
    !Generate the index for inequivalent sites 
    mag_basis=0
    mag_basis(1)=1
    do i=1,nsite
       r1=xyz(i,:)
       if(mag_basis(i).eq.0) then
          write(*,*) 'error' 
          stop
       endif
       jloop:do j=1,3
          r2=r1+d(j,:)
          call findsite(r2,id)
          !write(*,*) i,id
          if(id.eq.-1) cycle jloop
          select case (j)
          case(1,3)
             mag_basis(id)=mag_basis(i)+1
             if(mag_basis(id).gt.3) mag_basis(id)=1
          case(2)
             mag_basis(id)=mag_basis(i)-1
             if(mag_basis(id).lt.1) mag_basis(id)=3
          end select
       enddo jloop
    enddo
  end subroutine identify_magnetic_basis

         
  subroutine findsite(r,id)
    implicit none
    real(kind=double),dimension(:),intent(in)::r
    integer,intent(out)::id
    integer::i,j,f
    real(kind=double)::dis
    f=-1
    do i=1,nsite
       dis=(r(1)-xyz(i,1))**2 + (r(2)-xyz(i,2))**2 
       if(dis.lt.1d-5) then
         id=i; f=0
         exit
       endif
    enddo 
    if(f.ne.0) id=-1
  end subroutine findsite
    
end module magnetic_basis
  

