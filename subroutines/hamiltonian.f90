!  Module: hamiltonian
!
!  This modules generates and optionally diagonalizes the single particle Hamiltonian 
!  on any lattice (defined in NN table)  with a given onsite potential disorder
!
!  Author: Oinam Nganba Meetei
!  Date:   02/08/11
!  
!  Modifications: 08/15/12
!  Twisted BC is added.


module hamiltonian

  use mtprng

  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer,dimension(:,:),allocatable,save::NN
  integer,save::nsite
  type(mtprng_state),save::state
  integer,save::init_index=-1

  public:: init_ham,generate_disorder,generate_hamiltonian,diagonalize,state

  !Define common subroutine interface names
  interface diagonalize
     module procedure diagonalize_complex
     module procedure diagonalize_real
  end interface diagonalize
  interface generate_hamiltonian
     module procedure generate_real_hamiltonian
     module procedure generate_complex_hamiltonian
  end interface generate_hamiltonian




contains
  subroutine init_ham()
    implicit none
    integer::i,j,INFO
    open(unit=11,file='NN.out')
    read(11,*) nsite
    allocate(NN(nsite,20),stat=INFO)
    if(INFO.ne.0) then
       write(*,*) 'Allocation of neighbors array failed'
       stop
    endif
    do i=1,nsite
       read(11,*) NN(i,1)
       read(11,*)(NN(i,j),j=2,NN(i,1)+1)
    enddo
    close(11)
    write(*,*) 'nsite',nsite
    init_index=0
    write(*,*) 'hamiltonian initialized'
  end subroutine init_ham
  


  subroutine generate_disorder(dV,seed,V)
    !Generate disorder uniformly from the range [-V/2,V/2]
    implicit none
    real(kind=double),intent(in)::dV
    integer,intent(in)::seed
    real(kind=double),dimension(:),intent(out)::V
    integer::i
    if(init_index.ne.0) then
       write(*,*) 'hamiltonian not initialized yet'
       stop
    endif
    if(size(V).ne. nsite) then
       write(*,*) 'Dimension of V doesnt match system size'
       stop
    endif
    call mtprng_init(seed,state)
    do i=1,nsite
       V(i)=dV*(mtprng_rand_real3(state)-0.5)
    enddo
    V=V-sum(V)/nsite
    if(abs(sum(V)).gt.1d-5) then 
       write(*,*) 'sum V',sum(V)
       stop
    endif
    write(*,*) 'disorder generated'
    !write(*,*) V
  end subroutine generate_disorder


  subroutine generate_real_hamiltonian(t,mu,V,hmat)
    implicit none
    real(kind=double),intent(in)::t,mu
    real(kind=double),dimension(:),intent(in)::V
    real(kind=double),dimension(:,:),intent(out)::hmat
    integer::i,j
    if(init_index.ne.0 .or. size(V).ne.nsite .or. size(hmat,1).ne.nsite) then
       write(*,*) 'Error in generating Hamiltonian'
       stop
    endif
    !diagonal part
    hmat=0.0
    do i=1,nsite
       hmat(i,i)=V(i)-mu
    enddo
    !NN hopping
    do i=1,nsite
       do j=1,NN(i,1)
          hmat(i,NN(i,j+1))=-t
       enddo
    enddo
    write(*,*) 'Hamiltonian generated'
  end subroutine generate_real_hamiltonian


  subroutine generate_complex_hamiltonian(t,mu,V,hmat,theta,phi)
    implicit none
    real(kind=double),intent(in)::t,mu
    real(kind=double),optional::theta,phi
    real(kind=double),dimension(:),intent(in)::V
    complex(kind=double),dimension(:,:),intent(out)::hmat
    integer::i,j,k,Nx,Ny,Nz,m1,m2,m3
    real(kind=double),dimension(3,3)::latvec
    real(kind=double),dimension(3)::r
    real(kind=double),dimension(nsite)::x,y,z
    real(kind=double)::r_NN,dis
    real(kind=double)::th,ph
    if(init_index.ne.0 .or. size(V).ne.nsite .or. size(hmat,1).ne.nsite) then
       write(*,*) 'Error in generating Hamiltonian'
       stop
    endif
    !Read lattice parameters (size and lattice vectors)
    open(unit=11,file='square_lat_param.inp',status='OLD')
    read(11,*) Nx,Ny,Nz
    do i=1,3
       read(11,*)(latvec(i,j),j=1,3)
    enddo
    close(11)
    !Read co-ordinates
    open(unit=11,file='xyz.out',status='OLD')
    do i=1,nsite
       read(11,*) x(i),y(i),z(i)
    enddo   
    close(11) 
    if(.not.present(theta) .or. .not.present(phi)) then
       th=0d0; ph=0d0
    else
       th=theta;ph=phi
    endif
    !diagonal part
    hmat=0.0
    do i=1,nsite
       hmat(i,i)=V(i)-mu
    enddo
    !NN hopping
    r_NN=1.d0
    do i=1,nsite
       do k=2,NN(i,1)+1
          j=NN(i,k)
          hmat(i,j)=-t
          !check if the bond wraps around the lattice and add phase twist
          do m1=-1,1
          do m2=-1,1
             if(m1.eq.0 .and. m2.eq.0) cycle
             r(1) = x(i)-x(j)-m1*Nx*latvec(1,1)-m2*Ny*latvec(2,1)!-m3*Nz*basis(3,1)
             r(2) = y(i)-y(j)-m1*Nx*latvec(1,2)-m2*Ny*latvec(2,2)!-m3*Nz*basis(3,2)
             r(3) = z(i)-z(j)-m1*Nx*latvec(1,3)-m2*Ny*latvec(2,3)!-m3*Nz*basis(3,3)
             dis = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
             if(abs(dis-r_NN).le.1d-2) then
                select case (abs(m1*m2))
                case (0)
                   !write(*,*) x(i)-x(j), y(i)-y(j)
                   if(m2.eq.0) then
                      if(x(i).lt.x(j)) then
                         hmat(i,j)=-t*exp((0d0,-1.d0)*th)
                      else
                         hmat(i,j)=-t*exp((0d0,1.d0)*th)
                      endif
                   else
                      if(y(i).lt.y(j)) then
                         hmat(i,j)=-t*exp((0d0,-1.d0)*ph)
                      else
                         hmat(i,j)=-t*exp((0d0,1.d0)*ph)
                      endif
                      
                   endif
                   !write(*,*) hmat(i,j)
                case (1)
                   !write(*,*) x(i)-x(j), y(i)-y(j)
                   if(x(i).lt.x(j) .and. y(i).lt.y(j)) then
                      hmat(i,j)=-t*exp((0d0,1.d0)*(-th-ph))
                   elseif(x(i).lt.x(j) .and. y(i).gt.y(j)) then
                      hmat(i,j)=-t*exp((0d0,1.d0)*(-th+ph))
                   elseif(x(i).gt.x(j) .and. y(i).lt.y(j)) then
                      hmat(i,j)=-t*exp((0d0,1.d0)*(th-ph)) 
                   elseif(x(i).gt.x(j) .and. y(i).gt.y(j)) then
                      hmat(i,j)=-t*exp((0d0,1.d0)*(th+ph))
                   endif
                   !write(*,*) hmat(i,j)  
                end select
             endif
          enddo
          enddo
       enddo
    enddo
    !write(*,*) 'Hamiltonian generated'
  end subroutine generate_complex_hamiltonian
  


  subroutine diagonalize_real(ham,eval)
    implicit none
    real(kind=double),dimension(:,:),intent(inout)::ham
    real(kind=double),dimension(:),intent(out)::eval
    integer::LWORK,ndim
    real(kind=double),dimension(:),allocatable::WORK
    integer::INFO
    ndim=size(ham,1)
    LWORK=3*ndim
    allocate(WORK(LWORK),stat=INFO)
    if(INFO.ne.0) then
       write(*,*) 'allocation of diagonalization arrays failed'
       stop
    endif
    call DSYEV('V', 'U', ndim, ham, ndim, eval, WORK, LWORK, INFO )
    if(INFO.ne.0) then
       write(*,*) 'diagonalization failed'
       stop
    endif
    deallocate(WORK)
    !write(*,*) 'Hamiltonian diagonalized'
  end subroutine diagonalize_real

  

  subroutine diagonalize_complex(ham,eval)
    implicit none
    complex(kind=double),dimension(:,:),intent(inout)::ham
    real(kind=double),dimension(:),intent(out)::eval
    integer::LWORK,ndim
    complex(kind=double),dimension(:),allocatable::WORK
    real(kind=double),dimension(:),allocatable::RWORK
    integer::INFO
    ndim=size(ham,1)
    LWORK=3*ndim
    allocate(WORK(LWORK),RWORK(3*ndim-2),stat=INFO)
    if(INFO.ne.0) then
       write(*,*) 'allocation of diagonalization arrays failed'
       stop
    endif
    call ZHEEV('V', 'U', ndim, ham, ndim, eval, WORK, LWORK, RWORK, INFO )
    if(INFO.ne.0) then
       write(*,*) 'diagonalization failed'
       stop
    endif
    deallocate(WORK,RWORK)
    !write(*,*) 'Hamiltonian diagonalized'
  end subroutine diagonalize_complex


end module hamiltonian
    
