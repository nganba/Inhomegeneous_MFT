module sort_util
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  public::sort,sortby

  interface sort
     module procedure sort_vector
     module procedure sort_array
  end interface sort
  
contains
  
  subroutine sort_vector(x,n,o)
    implicit none 
    real(kind=double),dimension(:),intent(inout)::x
    integer,intent(in)::n
    character(1),intent(in),optional::o
    integer::i
    character(1)::choice
    real(kind=double),dimension(n)::tmpvec
    tmpvec=x(1:n)
    choice='a'
    if(present(o)) choice=o
    if(choice=='a') then
      call HPSORT(n,tmpvec)
      x(1:n)=tmpvec
    elseif(choice=='d') then
      call HPSORT(n,tmpvec)
      do i=0,n-1
        x(i+1)=tmpvec(n-i)
      enddo
    else
      write(*,*) 'Provide a or d as options'
      stop
    endif   
  end subroutine sort_vector


  subroutine sort_array(x,m,n,o)
    implicit none
    real(kind=double),dimension(:,:),intent(inout)::x
    integer,intent(in)::m,n
    character(1),intent(in),optional::o
    integer::i,j
    character(1)::choice
    real(kind=double),dimension(m,n)::tmparray
    real(kind=double),dimension(m)::tmpvec
    tmparray=x(1:m,1:n)
    if(present(o)) then
      choice=o
    else
      choice='a'
    endif
    if(choice=='a') then
      do i=1,n
         tmpvec=tmparray(:,i)
         call HPSORT(m,tmpvec) 
         x(1:m,i)=tmpvec
      enddo
    elseif(choice=='d') then
      do i=1,n
         tmpvec=tmparray(:,i)
         call HPSORT(m,tmpvec)
         do j=0,m-1
            x(j+1,i)=tmpvec(m-j)
         enddo
      enddo 
    else
      write(*,*) 'Provide a or d as options'
      stop
    endif
  end subroutine sort_array


  subroutine sortby(x,m,n,col,o)
    implicit none
    real(kind=double),dimension(:,:),intent(inout)::x
    integer,intent(in)::m,n,col
    character(1),intent(in),optional::o
    real(kind=double),dimension(m,n)::tmparray
    real(kind=double),dimension(m)::tmpvec
    integer,dimension(m)::indx
    character(1)::choice
    integer::i,j
    tmparray=x(1:m,1:n)
    choice='a'
    do i=1,m
       indx(i)=i
    enddo
    if(present(o)) choice=o
    if(choice=='a') then
       tmpvec=tmparray(:,col)
       call HPSORT2(m,tmpvec,indx)
       do i=1,m
          x(i,1:n)=tmparray(indx(i),:)
       enddo
    elseif(choice=='d') then
       tmpvec=tmparray(:,col)
       call HPSORT2(m,tmpvec,indx)
       do i=0,m-1
          x(i+1,1:n)=tmparray(indx(m-i),:)
       enddo
    else
       write(*,*) 'Provide a or d as options'
       stop
    endif
  end subroutine sortby 
       

  
           
  subroutine HPSORT(N,RA)
    !*****************************************************
    !*  Sorts an array RA of length N in ascending order *
    !*                by the Heapsort method             *
    !* ------------------------------------------------- *
    !* INPUTS:                                           *
    !*          N     size of table RA                   *
    !*          RA    table to be sorted                 *
    !* OUTPUT:                                           *
    !*          RA    table sorted in ascending order    *
    !*                                                   *
    !* NOTE: The Heapsort method is a N Log2 N routine,  *
    !*       and can be used for very large arrays.      *
    !*****************************************************  
    integer,intent(IN)::N
    real(kind=double),dimension(N),intent(inout)::RA
    integer::L,IR,I,J
    real(kind=double)::RRA
    L=N/2+1
    IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
10  continue
    if(L > 1)then
      L=L-1
      RRA=RA(L)
    else
      RRA=RA(IR)
      RA(IR)=RA(1)
      IR=IR-1
      if(IR.eq.1)then
        RA(1)=RRA
        return
      end if
    end if
    I=L
    J=L+L
20  if(J.le.IR)then
      if(J < IR)then
        if(RA(J) < RA(J+1))  J=J+1
      end if
      if(RRA < RA(J))then
        RA(I)=RA(J)
        I=J; J=J+J
      else
        J=IR+1
      end if
      goto 20
    end if
    RA(I)=RRA
    goto 10
  end subroutine HPSORT


  subroutine HPSORT2(N,RA,INDX)
    !*****************************************************
    !*  Sorts an array RA of length N in ascending order *
    !*  by the Heapsort method while tracking change in  *
    !*  index                                            *
    !* ------------------------------------------------- *
    !* INPUTS:                                           *
    !*          N       size of table RA                 *
    !*          RA      vector to be sorted              *
    !*          INDX    index of original vector         *
    !* OUTPUT:                                           *
    !*          RA      vector sorted in ascending order *
    !*          INDX    index of sorted vector           *
    !*                                                   *
    !* NOTE: The Heapsort method is a N Log2 N routine,  *
    !*       and can be used for very large arrays.      *
    !*****************************************************  
    integer,intent(IN)::N
    real(kind=double),dimension(N),intent(inout)::RA
    real(kind=double)::RRA
    integer,dimension(N)::INDX
    integer::L,IR,I,J,NINDX
    L=N/2+1
    IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
10  continue
    if(L > 1)then
      L=L-1
      RRA=RA(L)
      NINDX=INDX(L)
    else
      RRA=RA(IR)
      NINDX=INDX(IR)
      RA(IR)=RA(1)
      INDX(IR)=INDX(1)
      IR=IR-1
      if(IR.eq.1)then
        RA(1)=RRA
        INDX(1)=NINDX
        return
      end if
    end if
    I=L
    J=L+L
20  if(J.le.IR)then
      if(J < IR)then
        if(RA(J) < RA(J+1))  J=J+1
      end if
      if(RRA < RA(J))then
        RA(I)=RA(J)
        INDX(I)=INDX(J)
        I=J; J=J+J
      else
        J=IR+1
      end if
      goto 20
    end if
    RA(I)=RRA
    INDX(I)=NINDX
    goto 10
  end subroutine HPSORT2

end module sort_util
