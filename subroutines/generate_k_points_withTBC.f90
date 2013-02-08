! Program: generate_k_points_with_TBC.f90
! 
! From the given primitive vectors in real space, this code generates the 
! k-points inside the first Brillouin zone with a twisted boundary condition
!
! Author: Oinam Nganba Meetei
! Date:  08/16/12
!
!**************************************************************************

module k_points

  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.d0)
  integer,save :: Lx,Ly,Lz,n,init=-1
  real(kind=double),save ::a,b
  real(kind=double),dimension(3),save::a1,a2,a3,b1,b2,b3
  real(kind=double),dimension(20,3),save::NN
  
  public:: init_kvec,generate_kpoints
  
contains

  subroutine init_kvec(basis)
    implicit none
    real(kind=double),dimension(:,:),intent(in),optional::basis
    integer::i,j,k
    real(kind=double)::kx,ky,kz,dis
    
    !Read primitive vectors of real space lattice
    open(unit=22,file='triangular_lat_param.inp')
    read(22,*) Lx,Ly,Lz
    if(present(basis)) then
       do j=1,3
          a1(j)=basis(1,j)
          a2(j)=basis(2,j)
          a3(j)=basis(3,j)
       enddo
    else 
       read(22,*)(a1(j),j=1,3)
       read(22,*)(a2(j),j=1,3)
       read(22,*)(a3(j),j=1,3)
       init=0
    endif
    close(22)
    write(*,*) 'Primitive vectors'
    write(*,*) (a1(i),i=1,3)
    write(*,*) (a2(i),i=1,3)
    write(*,*) (a3(i),i=1,3)

    ! Calculate reciprocal lattice vectors
    b1(1)=a2(2)*a3(3)-a2(3)*a3(2)
    b1(2)=a2(3)*a3(1)-a2(1)*a3(3)
    b1(3)=a2(1)*a3(2)-a2(2)*a3(1)
    b1=((2*pi)/(dot_product(a1,b1)))*b1

    b2(1)=a3(2)*a1(3)-a3(3)*a1(2)
    b2(2)=a3(3)*a1(1)-a3(1)*a1(3)
    b2(3)=a3(1)*a1(2)-a3(2)*a1(1)
    b2=((2*pi)/(dot_product(a2,b2)))*b2

    b3(1)=a1(2)*a2(3)-a1(3)*a2(2)
    b3(2)=a1(3)*a2(1)-a1(1)*a2(3)
    b3(3)=a1(1)*a2(2)-a1(2)*a2(1)
    b3=((2*pi)/(dot_product(a3,b3)))*b3

    b=sqrt(b1(1)**2 + b1(2)**2 + b1(3)**2)

    write(*,*)
    write(*,*) 'Reciprocal lattice vectors'
    write(*,*) (b1(i),i=1,3)
    write(*,*) (b2(i),i=1,3)
    write(*,*) (b3(i),i=1,3)
 
    !Find nearest neighbors in reciprocal lattice
    !Assumes that only on length 'a' defines NN sites uniquely
    NN=0d0
    n=0
    do  i=-1,1
       do j=-1,1
          do k=-1,1
             kx=i*b1(1) + j*b2(1) + k*b3(1)
             ky=i*b1(2) + j*b2(2) + k*b3(2)
             kz=i*b1(3) + j*b2(3) + k*b3(3)
             dis=sqrt(kx**2 + ky**2 + kz**2)
             if(abs(dis-b).lt.1d-2) then
                n=n+1
                NN(n,1)=kx
                NN(n,2)=ky
                NN(n,3)=kz
             endif
          enddo
       enddo
    enddo

  end subroutine init_kvec


  subroutine generate_kpoints(theta,phi,bz,nkpts)
    implicit none
    integer,intent(out)::nkpts
    real(kind=double),intent(in)::theta,phi
    real(kind=double),dimension(:,:),intent(out)::bz
    integer::i,j,k,iter
    real(kind=double)::kx,ky,kz,dis,dis2
    real(kind=double),dimension(3)::ap1,ap2,ap3,bp1,bp2,bp3
 
    !Read original lattice parameters in case we are generating BZ for a
    !different lattice parameter
    if(init.ne.0) then
       open(unit=22,file='triangular_lat_param.inp')
       read(22,*) Lx,Ly,Lz
       read(22,*)(ap1(j),j=1,3)
       read(22,*)(ap2(j),j=1,3)
       read(22,*)(ap3(j),j=1,3)

       bp1(1)=ap2(2)*ap3(3)-ap2(3)*ap3(2)
       bp1(2)=ap2(3)*ap3(1)-ap2(1)*ap3(3)
       bp1(3)=ap2(1)*ap3(2)-ap2(2)*ap3(1)
       bp1=((2*pi)/(dot_product(ap1,bp1)))*bp1

       bp2(1)=ap3(2)*ap1(3)-ap3(3)*ap1(2)
       bp2(2)=ap3(3)*ap1(1)-ap3(1)*ap1(3)
       bp2(3)=ap3(1)*ap1(2)-ap3(2)*ap1(1)
       bp2=((2*pi)/(dot_product(ap2,bp2)))*bp2

       bp3(1)=ap1(2)*ap2(3)-ap1(3)*ap2(2)
       bp3(2)=ap1(3)*ap2(1)-ap1(1)*ap2(3)
       bp3(3)=ap1(1)*ap2(2)-ap1(2)*ap2(1)
       bp3=((2*pi)/(dot_product(ap3,bp3)))*bp3
    else
       bp1=b1; bp2=b2; bp3=b3
    endif
   
    nkpts=0
    bz=0d0
    do i=-Lx,Lx
       do j=-Ly,Ly
          inner:do k=-Lz,Lz
             kx=(i + (theta/(2*pi)))*bp1(1)/Lx + (j + (phi/(2*pi)))*bp2(1)/Ly + k*bp3(1)/Lz
             ky=(i + (theta/(2*pi)))*bp1(2)/Lx + (j + (phi/(2*pi)))*bp2(2)/Ly + k*bp3(2)/Lz
             kz=(i + (theta/(2*pi)))*bp1(3)/Lx + (j + (phi/(2*pi)))*bp2(3)/Ly + k*bp3(3)/Lz
             dis=sqrt(kx**2 + ky**2 + kz**2)
             !Check if k belongs to 1st BZ
             if(abs(kz).gt.1d-4) cycle inner ! use this for 2D lattice
             do iter=1,n
                dis2=sqrt((kx-NN(iter,1))**2 + (ky-NN(iter,2))**2 + (kz-NN(iter,3))**2)
                if((dis-dis2).gt.1d-2) cycle inner
             enddo
             nkpts=nkpts+1
             bz(nkpts,1)=kx; bz(nkpts,2)=ky; bz(nkpts,3)=kz
          enddo inner
       enddo
    enddo

  end subroutine generate_kpoints
  
end module k_points
  
  
