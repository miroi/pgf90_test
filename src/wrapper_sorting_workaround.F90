  module wrapper_sorting

    implicit none

    public shuffle_indices
    public match_array
    public onebody_for_blas
    public sorting_for_blas
    public decide_sorting_routine 
    public sorted_dimension 

    private

    integer      :: sortType
    logical      :: do_direct, do_exchange
    logical      :: t_direct, t_exchange
    character(20):: arraytype
    character(2),pointer :: array_original(:)
    character(2)         :: exchange(4)
    integer              :: reorder_left(4) 
    type sorting_list
      integer :: num 
      integer :: trans
!     contains
      procedure(),pointer, nopass :: unrestricted !, restricted 
      procedure(),pointer, nopass :: restricted_braket, restricted_bra, restricted_ket 
    endtype sorting_list 

    type indices
      integer, dimension(:,:),pointer :: tuple
      integer, dimension(:),pointer   :: tuple_combInd, orbital_type, full_array 
    endtype indices

    type, extends(indices) :: middle_indices
       type(indices) :: triangular
    endtype

   type, extends(middle_indices) :: more_indices
      type(middle_indices) :: pqr,qrs,totsym
   end type

    type bit_processing
       integer :: bra,ket,braKet,I,J,K,L,ijk,jkl,ijkl
       integer :: I_reorder,J_reorder,K_reorder,L_reorder,ijkl_reorder
       integer :: bra_reorder,ket_reorder,ijk_reorder,jkl_reorder
    endtype bit_processing 

    type presorting
       real*8, allocatable :: presorted(:) 
       real*8, pointer :: integral(:)
    endtype
   external SRT1SS4
   external SRT1L1
   external SRT1C1
   external SRT16
   external SRT1R1
   external SRT1TT4
   external SRT1TT5
   external SRT1ST4
   external SRT36 
   external SRT22
   external SRT32
   external SRT26
   external SRT20D
   external SRT1T3
   external SRT1T2
   external SRT1S2
   external SRT19
   external SRT1S3
   external SRT7
   external SRT6
   external SRT9
   external SRT1TS4

   external getvoov
   external getoovv
   external getooov
   external getovoo
   external getoovo
   external getvvvo_incore


    contains
    
   subroutine shuffle_indices(reorder_array,no_match,save_index,trans)

! ------- variables --------
    integer,intent(inout)        :: reorder_array(4)
    logical,intent(inout)        :: no_match
    character(1),intent(inout)   :: trans
    integer                      :: shifted_array(4)
    integer,allocatable          :: P(:,:)
    integer                      :: temp,i,j  
    integer, intent(out)         :: save_index
    type(sorting_list)           :: a
! --------------------------

  save_index = 0

  select case(sortType)
    
    case(13)

     shifted_array(1:3) =  reorder_array(2:4)

     temp = reorder_array(1)

    allocate (P(6,3))

    call permutate(shifted_array(1:3),P)

    do i = 1, size(P,1)
      a%num = toDecimal((/temp,P(i,:)/))
      a%trans = toDecimal((/temp,P(i,:)/),.true.)    

      call match_array (a%num, no_match, save_index, trans)  
      if (no_match) call match_array_transpose (a%trans, no_match, save_index, trans)  

      if (.not.no_match) then
        reorder_array(2:4)=P(i,:)
        exit
      endif

    enddo
    deallocate(P)

     case(31) 

     shifted_array(1:3) = reorder_array(1:3)

     temp = reorder_array(4)

    allocate (P(6,3))

    call permutate(shifted_array(1:3),P)

    do i = 1, size(P,1)
!    write(*,*)P(i,:)
    a%num = toDecimal((/P(i,:),temp/))
    a%trans = toDecimal((/P(i,:),temp/),.true.)    

    call match_array (a%num, no_match, save_index, trans)  
    if (no_match) call match_array_transpose (a%trans, no_match, save_index, trans)  
    if (.not.no_match) then
    reorder_array(1:3)=P(i,:)
    exit
    endif
    enddo
    deallocate(P)

   case (22) 

     shifted_array(1:4) = reorder_array(1:4) 
     i = 1
   do
     temp = reorder_array(i)
     reorder_array(i) = reorder_array(i+1)
     reorder_array(i+1) = temp

    a%num = toDecimal(reorder_array)
    a%trans = toDecimal(reorder_array,.true.)    

    call match_array (a%num,no_match,save_index,trans)  
    if (no_match) call match_array_transpose (a%trans,no_match,save_index,trans) 
    if ((.not.no_match).or.(i>3)) exit
    i = i+2
    reorder_array(1:4) = shifted_array(1:4)
   enddo
  end select

   end subroutine

  subroutine match_array (permute,no_match,save_index,trans) 

!-----------descripton---
! choose the sorting type depending on whether the array is square or tringular and pattern of sorting (e.g, (2,2);(1,3))
!------------------------

!-----------calling variable-------------
   integer, intent(in)        :: permute 
   logical, intent(inout)     :: no_match
   integer, intent(out)       :: save_index
   character(1), intent(inout):: trans
!--------------------------------
!---------local variable---------
   integer, allocatable       :: A(:)
   integer                    :: i,dim_A,upbound_A
!--------------------------------

  save_index = 0

  select case(arraytype)

  case('unrestricted')

  select case(sortType)

  case(22) 

  allocate(A(5))
  
  A = (/1234,1324,1432,2134,1243/)

  case(31)

  allocate(A(5))

  A = (/1234,1243,1342,3421,4123/)

  case(13)

  allocate(A(1))

  A = (/1234/)

  end select

  case('restricted_braket')

  select case(sortType) 

  case(22)
  allocate(A(3))
  A = (/1234,1324,1342/)

  case(31)
 
  allocate(A(3))

  A= (/1234,4123,1342/)
  
  case(13) 
  allocate(A(1))

  A= (/1234/)

  end select
 
  case('restricted_bra')

  select case(sortType) 

  case(22)
  allocate(A(3))
  A = (/1234,1432,1324/)

  case(31)
  allocate(A(4))
  A = (/1234,1342,1243,4123/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select
 
  case('restricted_ket')

  select case(sortType)

  case(22)

  allocate(A(3))
  A = (/1234,1324,1432/)

  case(31)

  allocate(A(4))
  A = (/1234,4123,1342,3421/)


  case(13)
  allocate(A(1))
  A = (/1234/)

  end select

  end select

  dim_A = size(A)
  upbound_A = dim_A

   i = 1

   do 
   if (permute == A(i)) then
   save_index = i
   no_match = .FALSE.
   exit
   endif
     i = i+1
   if (i > upbound_A) exit
   enddo 

  deallocate(A)

  end subroutine

  subroutine match_array_transpose (permute_transpose, no_match, save_index, trans) 

!-----------variable-------------
   integer,intent(in) :: permute_transpose
   logical,intent(inout):: no_match
   integer,intent(out):: save_index
   character(1),intent(inout) :: trans
!--------------------------------
!---------local variable---------
  integer, allocatable :: A(:)
  integer :: i, dim_A, upbound_A
  integer :: sort_type
!--------------------------------

  save_index = 0 

  if (sortType == 31) sort_type = 13
  if (sortType == 13) sort_type = 31
  if (sortType == 22) sort_type = 22

 select case(arraytype)

  case('unrestricted')

  select case(sort_type)

  case(22) 

  allocate(A(5))

  A = (/1234,1324,1432,2134,1243/)

  case(31)

  allocate(A(5))

  A = (/1234,1243,1342,3421,4123/)

  case(13)

  allocate(A(1))

  A = (/1234/)

  end select

  case('restricted_braket')

  select case(sort_type) 

  case(22)
  allocate(A(3))
  A = (/1234,1324,1342/)

  case(31)
 
  allocate(A(3))

  A= (/1234,4123,1342/)
  
  case(13) 
  allocate(A(1))

  A= (/1234/)

  end select
 
  case('restricted_bra')

  select case(sort_type) 

  case(22)
  allocate(A(3))
  A = (/1234,1432,1324/)

  case(31)
  allocate(A(4))
  A = (/1234,1342,1243,4123/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select
 
  case('restricted_ket')

  select case(sort_type)

  case(22)

  allocate(A(3))
  A = (/1234,1324,1432/)

  case(31)

  allocate(A(4))
  A = (/1234,4123,1342,3421/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select

  end select

  dim_A = size(A)
  upbound_A = dim_A

   i = 1
   do 
    if (permute_transpose == A(i)) then 
      save_index = i
      no_match = .FALSE.
!      trans ='C'  
       trans ='T'  
      exit
    endif
     i = i+1
    if (i > upbound_A) exit
   enddo 

  deallocate(A)

  end subroutine

  subroutine decide_sorting_routine(sort_type,opArray_left,opArray_right,opArray_target,if_trans, &
                                &    save_index,asym_factor,do_shuffle,sorted_opArray,inverse)

!-------------------------------------
!    choose the sorting routine. 
!-------------------------------------

!--------------calling variables---------------
    character(2), intent(in),target        :: opArray_left(:)
    character(2), intent(in)               :: opArray_target(:),opArray_right(:)
    integer,      intent(in)               :: sort_type
    logical,      intent(in), optional     :: inverse  
    character(1), intent(out)              :: if_trans
    integer,      intent(out)              :: save_index 
    integer,      intent(out)              :: asym_factor
    logical,      intent(in)               :: do_shuffle
    character(2), intent(out),optional     :: sorted_opArray(size(opArray_left))
!----------------------------------------------

!-------local variables--------
    character(2) :: array_exchange(size(opArray_left))
    integer      :: reorder_v(size(opArray_left))
    integer      :: i,j,k,l 
    logical      :: no_match  
    character(2) :: p,q,r,s

    type(sorting_list)     :: a
!------------------------------

     asym_factor = 1

     sortType = sort_type       !very dirty trick. will soon amend it.

     array_original => opArray_left 

     arraytype = 'unrestricted'

     p = opArray_left(1)
     q = opArray_left(2)
     r = opArray_left(3)
     s = opArray_left(4)

    if (p(1:1) == q(1:1)) then
      if (r(1:1) == s(1:1)) then
       arraytype = 'restricted_braket'
      else
       arraytype = 'restricted_bra'
      endif 
      elseif (r(1:1)==s(1:1)) then
       arraytype = 'restricted_ket'
    endif

    do_direct = .true.
    do_exchange = .false.
    no_match = .true.

   if (present(inverse)) then 

    call align_array(array1=opArray_left,array2=opArray_target,reorder_array=reorder_v)

!findloc feature is not available for intel compilers. so I can't use it now.

    a%num = toDecimal(reorder_v)    
    a%trans = toDecimal(reorder_v,.true.)    

    call match_array (a%num,no_match,save_index,if_trans)

    if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

    endif

    if ((no_match).and.(do_direct)) then

    call align_array(opArray_left,opArray_target,opArray_right,reorder_v)

!findloc feature is not available for intel compilers. so I can't use it now.

    a%trans = toDecimal(reorder_v,.true.)    
    a%num   = toDecimal(reorder_v)    

    call match_array (a%num,no_match,save_index,if_trans)

    if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  


!   if (arraytype=='unrestricted') then 

    if ((no_match).and.(do_shuffle)) call shuffle_indices (reorder_v,no_match,save_index,if_trans)

!   endif

     if (no_match) then
      write (*,*) "====no matching sorting subroutine available for direct mode of contraction.====="
      do_exchange = .true.
     endif

    endif 

   if (do_exchange) then

   !! In this mode I generate all the antisymmetrized components of an array on the fly where the interchange is possible between the same orbital-type (i.e, h or p).

   !! i have included the change in sign due to anti-symmetrization.

       write(*,*)'=====contraction will be done in exchange mode====='

       array_exchange = opArray_left

     if (arraytype == 'restricted_braket') then
        k = 1
        i = 1 

     do

       call swap_character(array_exchange,i,i+1)   

       call align_array(array_exchange,opArray_target,opArray_right,reorder_v)

       a%num   = toDecimal(reorder_v)
       a%trans = toDecimal(reorder_v,.true.)    

       call match_array (a%num,no_match,save_index,if_trans)
      
       if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

       if (.not.no_match) exit
       if (k>3)           exit
       i = i + size(opArray_left)/2 
       if (k>2) then
        array_exchange = opArray_left
        i = 3
       endif
       k = k+1
     enddo

     if(.not.no_match) asym_factor = asym_factor*(-1)**k

!      if ((no_match).and.(do_shuffle)) call shuffle_indices (reorder_v,no_match,save_index,if_trans)

   else

       if (arraytype == 'restricted_ket') call swap_character(array_exchange,3,4)   
       if (arraytype == 'restricted_bra') call swap_character(array_exchange,1,2)   

       call align_array(array_exchange,opArray_target,opArray_right,reorder_v)

       a%num   = toDecimal(reorder_v)
       a%trans = toDecimal(reorder_v,.true.)    

       call match_array (a%num,no_match,save_index,if_trans)
      
       call match_array_transpose (a%trans,no_match,save_index,if_trans)  

       if ((no_match).and.(do_shuffle)) call shuffle_indices (reorder_v,no_match,save_index,if_trans)

       if(.not.no_match) asym_factor = -1*asym_factor

    endif
   endif

    exchange = array_exchange

!  if ((no_match).and.(do_shuffle)) then

!  do_exchange = .false.

!  call align_array(opArray_left,opArray_target,opArray_right,reorder_v)

!  call shuffle_indices (reorder_v,no_match,save_index,if_trans)

!  endif 


    if (no_match) then
      call quit('No Sorting subroutine available. code stops here.')
    endif

    if (present(sorted_opArray)) then 
      if (.not.do_exchange) then
        sorted_opArray = get_operator_array(opArray_left,toDecimal(reorder_v))
      else
        sorted_opArray = get_operator_array(array_exchange,toDecimal(reorder_v))
      endif
    endif

   reorder_left = reorder_v 

  if ((if_trans=='T').and.(sortType==31)) then

   reorder_left(1)=reorder_v(4)
   reorder_left(2:4)=reorder_v(1:3)

  endif 

  if ((if_trans=='T').and.(sortType==13)) then

   reorder_left(4)=reorder_v(1)
   reorder_left(1:3)=reorder_v(2:4)

  endif 

  if ((if_trans=='T').and.(sortType==22)) then

   reorder_left(1:2)=reorder_v(3:4)
   reorder_left(3:4)=reorder_v(1:2)

  endif 

  end subroutine

  subroutine sorting_for_blas(save_index,do_inverse,if_trans,sorted,row,column,presorted,if_TotSym) 

!---------------description-------
!        carry-out the sorting.
!---------------------------------

!--------------calling variables---------------
    real(8),intent(in),target,optional           :: presorted(:)
    logical,intent(in),optional           :: if_TotSym
    integer,intent(inout)                 :: save_index
    logical,intent(in)                    :: do_inverse  
    character(1),intent(inout)            :: if_trans
    integer,intent(out)                   :: row(:), column(:)
    real(8),intent(inout)                 :: sorted(:) !think a bit more carefully to replace it by dynamic allocation. not easy.
!----------------------------------------------

!-------local variables------------------------
!   character(2) :: array_exchange(4)
    integer      :: reorder_v(4)
    logical      :: do_TotSym
    integer      :: i,j,k,l 
    character(2) :: p,q,r,s
    integer      :: sort_type_dummy,shift,shift1,AllocateStatus
#include "symm.inc"
#include "param.inc"
#include "complex.inc"

    type(sorting_list)     :: a
    type(middle_indices)   :: c(0:30)  
    type(more_indices)     :: m(0:15)  
    type(bit_processing)   :: d 
    type(presorting),target  :: V 
    type(presorting)       :: B

    type(sorting_list) :: srt_22(0:5)
    type(sorting_list) :: srt_31(0:5)
    type(sorting_list) :: srt_13(0:5)
! (2,2)->(2,2) type of sorting routines         
       srt_22(1)%unrestricted => SRT1SS4   
       srt_22(2)%unrestricted => SRT16     
       srt_22(3)%unrestricted => SRT1L1    
       srt_22(4)%unrestricted => SRT1R1 
       srt_22(1)%restricted_braket =>  SRT1TT4
       srt_22(2)%restricted_braket =>  SRT1TT5
!      srt_22(1)%restricted_bra =  SRT1LS1
       srt_22(1)%restricted_bra => SRT36
       srt_22(2)%restricted_bra => SRT1TS4
       srt_22(1)%restricted_ket => SRT1ST4
       srt_22(2)%restricted_ket => SRT26
!       srt_22(3)%restricted_ket => SRT20D

! (2,2)->(3,1) type of sorting routines.
       srt_31(0)%unrestricted => SRT1S3
       srt_31(1)%unrestricted => SRT19
       srt_31(2)%unrestricted => SRT6
       srt_31(3)%unrestricted => SRT9
       srt_31(4)%unrestricted => SRT32
       srt_31(0)%restricted_ket => SRT1T3
       srt_31(1)%restricted_ket => SRT22
       srt_31(2)%restricted_ket => SRT6
       srt_31(3)%restricted_ket => SRT9
       srt_31(0)%restricted_braket => SRT1T3
       srt_31(1)%restricted_braket => SRT22
       srt_31(2)%restricted_braket => SRT7
       srt_31(0)%restricted_bra => SRT1S3
       srt_31(1)%restricted_bra => SRT7
       srt_31(2)%restricted_bra => SRT19
       srt_31(3)%restricted_bra => SRT32
    
!(2,2)->(1,3) type of sorting routines.
       srt_13(0)%unrestricted => SRT1S2
       srt_13(0)%restricted_bra => SRT1T2
       srt_13(0)%restricted_ket => SRT1S2
       srt_13(0)%restricted_braket => SRT1T2

!----------------------------------------------

     save_index = save_index - 1

      do_TotSym =.false.

     if (present(if_TotSym)) do_TotSym = if_TotSym

     shift     = 0
     shift1    = 0

     if (.not.present(presorted)) then
     call integral_fetching_incore(V) 

     B%integral => v%presorted(1:size(v%presorted)) 

     endif


     if (present(presorted)) then

      B%integral => presorted(1:size(presorted)) 

     endif

     sort_type_dummy = sortType

     if (if_trans == 'T') then
       if (sortType == 31) sort_type_dummy = 13
       if (sortType == 13) sort_type_dummy = 31
     endif

    if (sortType == 22) then

     if ((.not.do_TotSym).and.(.not.do_inverse)) shift=15  

     if (save_index > 0) then
      select case(arraytype)
        
         case('unrestricted')

       call input_index(c,d,m)

       if (do_inverse) then

       shift1 = 4

       if (if_trans=='N') then
        column =  c(d%ket+shift1)%tuple_combInd
        row    =  c(d%bra+shift1)%tuple_combInd
       else
        row    = c(d%ket+shift1)%tuple_combInd
        column = c(d%bra+shift1)%tuple_combInd
       endif
       else
       if (if_trans=='N') then
        row     = c(d%bra_reorder)%tuple_combInd
        column  = c(d%ket_reorder)%tuple_combInd     
       else 
        column  = c(d%bra_reorder)%tuple_combInd
        row     = c(d%ket_reorder)%tuple_combInd     
       endif
       endif

      if (save_index==3) then

       shift1 = 4
       call srt_22(save_index)%unrestricted (nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%I)%orbital_type, &
     &         c(d%J)%orbital_type,c(d%ket+shift1)%tuple_combInd,m(d%braket)%totsym%full_array,c(d%bra_reorder+shift1)%tuple, &
     &         B%integral,sorted) 

      elseif (save_index==4) then

       call srt_22(save_index)%unrestricted (nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%K)%orbital_type, &
     &         c(d%L)%orbital_type,c(d%ket+shift1)%tuple_combInd,m(d%braket)%totsym%full_array,c(d%ket_reorder+shift1)%tuple, &
     &         B%integral,sorted) 

      elseif (save_index==2) then

      call srt_22(save_index)%unrestricted (nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &     c(d%K)%orbital_type,c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braket+shift)%full_array, &
     &     c(d%bra_reorder)%tuple,c(d%ket_reorder)%tuple,B%integral,sorted)
    
      else

      call srt_22(save_index)%unrestricted (nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &     c(d%K)%orbital_type,c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braket+shift)%full_array, &
     &     c(d%bra_reorder)%tuple,c(d%ket_reorder)%tuple,B%integral,sorted)


       endif

         case('restricted_braket') 

      call input_index(c,d)

       if (do_inverse) then
       if (if_trans=='N') then
        row     = c(d%bra)%triangular%tuple_combInd
        column  = c(d%ket)%triangular%tuple_combInd     
       else 
        column  = c(d%bra)%triangular%tuple_combInd
        row     = c(d%ket)%triangular%tuple_combInd     
       endif
       else
       if (if_trans=='N') then
        row     = c(d%bra_reorder)%tuple_combInd
        column  = c(d%ket_reorder)%tuple_combInd     
       else 
        column  = c(d%bra_reorder)%tuple_combInd
        row     = c(d%ket_reorder)%tuple_combInd     
       endif
       endif

      
      if (save_index==2) then

      call srt_22(save_index)%restricted_braket(nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type,&
      &    c(d%K)%orbital_type,c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braket)%full_array, &
      &    c(d%bra_reorder)%tuple,c(d%ket_reorder)%tuple,B%integral,sorted)

      else

      call srt_22(save_index)%restricted_braket(nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type,&
      &    c(d%K)%orbital_type,c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braket+shift)%full_array, &
      &    c(d%bra_reorder)%tuple,c(d%ket_reorder)%tuple,B%integral,sorted)

      endif

         case('restricted_bra')

       call input_index(c,d)

       if (do_inverse) then
        shift1 = 4
       if (if_trans=='N') then
        row     = c(d%bra)%triangular%tuple_combInd
        column  = c(d%ket+shift1)%tuple_combInd     
       else 
        column  = c(d%bra)%triangular%tuple_combInd
        row     = c(d%ket+shift1)%tuple_combInd     
       endif
       else
       if (if_trans=='N') then
        row     = c(d%bra_reorder)%tuple_combInd
        column  = c(d%ket_reorder)%tuple_combInd     
       else 
        column  = c(d%bra_reorder)%tuple_combInd
        row     = c(d%ket_reorder)%tuple_combInd     
       endif
       endif

       call srt_22(save_index)%restricted_bra(nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
                  & c(d%K)%orbital_type,c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braket+shift)%full_array,  &
                  & c(d%bra_reorder)%tuple,c(d%ket_reorder)%tuple,B%integral,sorted)  

      case('restricted_ket')

      call input_index(c,d)

       if (do_inverse) then
        shift1 = 4
       if (if_trans=='N') then
        row     = c(d%bra+shift1)%tuple_combInd
        column  = c(d%ket)%triangular%tuple_combInd     
       else 
        column  = c(d%bra+shift1)%tuple_combInd
        row     = c(d%ket)%triangular%tuple_combInd     
       endif
       else
       if (if_trans=='N') then
        row     = c(d%bra_reorder)%tuple_combInd
        column  = c(d%ket_reorder)%tuple_combInd     
       else 
        column  = c(d%bra_reorder)%tuple_combInd
        row     = c(d%ket_reorder)%tuple_combInd     
       endif
       endif

      if (save_index==2) then

      call srt_22(save_index)%restricted_ket(nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type,  &
     &             c(d%K)%orbital_type,c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braket)%full_array,&
     &             c(d%bra_reorder)%tuple,c(d%ket_reorder)%tuple,B%integral,sorted)


      else

      call srt_22(save_index)%restricted_ket(nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type,  &
     &             c(d%K)%orbital_type,c(d%L)%orbital_type,c(d%bra_reorder)%tuple_combInd,c(d%braket+shift)%full_array,&
     &             c(d%bra_reorder)%tuple,c(d%ket_reorder)%tuple,B%integral,sorted)

      endif
       end select
     else
!!      write(*,*)'====sorting is not needed===='

     if (present(presorted)) then

     call dcopy(size(presorted),B%integral,1,sorted,1) 

     else

     call dcopy(size(V%presorted),B%integral,1,sorted,1) 

     endif

     call input_index(c,d)
       shift1 = 4
      select case(arraytype)
 
      case('restricted_braket')

       if (if_trans=='N') then
       row = c(d%bra)%triangular%tuple_combInd
       column = c(d%ket)%triangular%tuple_combInd     
       else 
       column = c(d%bra)%triangular%tuple_combInd
       row = c(d%ket)%triangular%tuple_combInd     
       endif

       case('restricted_bra')

       if (if_trans=='N') then
       row = c(d%bra)%triangular%tuple_combInd
       column = c(d%ket+shift1)%tuple_combInd     
       else 
       column = c(d%bra)%triangular%tuple_combInd
       row = c(d%ket+shift1)%tuple_combInd     
       endif

       case('restricted_ket')

       if (if_trans=='N') then
       row = c(d%bra+shift1)%tuple_combInd
       column = c(d%ket)%triangular%tuple_combInd     
       else 
       column = c(d%bra+shift1)%tuple_combInd
       row = c(d%ket)%triangular%tuple_combInd     
       endif

       case('unrestricted')

       if (if_trans=='N') then
       row    = c(d%bra+shift1)%tuple_combInd
       column = c(d%ket+shift1)%tuple_combInd     
       else 
       column = c(d%bra+shift1)%tuple_combInd
       row    = c(d%ket+shift1)%tuple_combInd     
       endif

       end select

     endif

     elseif (sort_type_dummy == 31) then
       select case(arraytype)
      
       case('unrestricted')

        call input_index(c,d,m)

       shift1 = 4

       if (save_index == 2) then

       if (if_trans=='N') then
        column = c(d%l_reorder)%orbital_type
        row = m(d%ijk_reorder)%qrs%tuple_combInd  
       else
        row = c(d%l_reorder)%orbital_type
        column = m(d%ijk_reorder)%qrs%tuple_combInd  
       endif
       
        call srt_31(save_index)%unrestricted(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%I)%orbital_type, &
           & c(d%J)%orbital_type,c(d%ket+shift1)%tuple_combInd,m(d%ijk_reorder)%qrs%tuple_combInd,&
           & array_tot(m(d%ijk_reorder)%qrs%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
           & m(d%ijk_reorder)%qrs%tuple,B%integral,sorted) 

        elseif (save_index == 3) then 

       if (if_trans=='N') then
        column = c(d%l_reorder)%orbital_type
        row = m(d%ijk_reorder)%pqr%tuple_combInd  
       else
        row = c(d%l_reorder)%orbital_type
        column = m(d%ijk_reorder)%pqr%tuple_combInd  
       endif
       
        call srt_31(save_index)%unrestricted(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%I)%orbital_type, &
           & c(d%J)%orbital_type,c(d%ket+shift1)%tuple_combInd,m(d%ijk_reorder)%pqr%tuple_combInd,&
           & array_tot(m(d%ijk_reorder)%pqr%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
           & m(d%ijk_reorder)%pqr%tuple,B%integral,sorted) 

        elseif (save_index == 4) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%qrs%tuple_combInd
       endif

     call srt_31(save_index)%unrestricted(nrep,multb,do_inverse,c(d%bra+shift1)%triangular%tuple_combInd,c(d%K)%orbital_type,&
           & c(d%L)%orbital_type,m(d%ijk_reorder)%qrs%tuple_combInd, &
           & array_tot(m(d%ijk_reorder)%qrs%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
           & m(d%ijk_reorder)%qrs%tuple,B%integral,sorted)

       else

       if (if_trans=='N') then
        column = c(d%l_reorder)%orbital_type
        row = m(d%ijk_reorder)%pqr%tuple_combInd  
       else
        row = c(d%l_reorder)%orbital_type
        column = m(d%ijk_reorder)%pqr%tuple_combInd  
       endif
       
        call srt_31(save_index)%unrestricted(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%K)%orbital_type, &
           & c(d%L)%orbital_type,m(d%ijk_reorder)%pqr%tuple_combInd,&
           & array_tot(m(d%ijk_reorder)%pqr%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
           & m(d%ijk_reorder)%pqr%tuple,B%integral,sorted) 

       endif

       if (do_inverse) then
       if (if_trans=='N') then
        column =  c(d%ket+shift1)%tuple_combInd
        row    = c(d%bra+shift1)%tuple_combInd
       else
        row = c(d%ket+shift1)%tuple_combInd
        column = c(d%bra+shift1)%tuple_combInd
       endif
       endif

         case('restricted_braket')

        call input_index(c,d,m)

       if (save_index == 2) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
       endif

      call srt_31(save_index)%restricted_braket(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%I)%orbital_type, &
        &  c(d%J)%orbital_type,m(d%ijk_reorder)%qrs%triangular%tuple_combInd, &
        &  array_tot(m(d%ijk_reorder)%qrs%triangular%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
        &  m(d%ijk_reorder)%qrs%triangular%tuple,B%integral,sorted,c(d%ket)%triangular%tuple_combInd)

       elseif (save_index == 1) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
       endif

     call srt_31(save_index)%restricted_braket(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type,&
           & c(d%L)%orbital_type,m(d%ijk_reorder)%qrs%triangular%tuple_combInd, &
           & array_tot(m(d%ijk_reorder)%qrs%triangular%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
           & m(d%ijk_reorder)%qrs%triangular%tuple,B%integral,sorted)

       else

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
       endif

     call srt_31(save_index)%restricted_braket(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type,&
              & c(d%L)%orbital_type,m(d%ijk_reorder)%pqr%triangular%tuple_combInd, &
              & array_tot(m(d%ijk_reorder)%pqr%triangular%tuple_combInd,c(d%l_reorder)%orbital_type,nrep),        &
              & m(d%ijk_reorder)%pqr%triangular%tuple,B%integral,sorted)

       endif

       if (do_inverse) then
       if (if_trans=='N') then
       row = c(d%bra)%triangular%tuple_combInd
       column = c(d%ket)%triangular%tuple_combInd     
       else 
       column = c(d%bra)%triangular%tuple_combInd
       row = c(d%ket)%triangular%tuple_combInd     
       endif
       endif

        case('restricted_bra')

        call input_index(c,d,m)

        shift1 = 4 
       
    if (save_index == 1) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%qrs%tuple_combInd
       endif

      call srt_31(save_index)%restricted_bra(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%I)%orbital_type, &
        &  c(d%J)%orbital_type,m(d%ijk_reorder)%qrs%tuple_combInd, &
        &  array_tot(m(d%ijk_reorder)%qrs%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
        &  m(d%ijk_reorder)%qrs%tuple,B%integral,sorted,c(d%ket+shift1)%tuple_combInd)

    elseif (save_index==3) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
       endif

     call srt_31(save_index)%restricted_bra(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type,&
           & c(d%L)%orbital_type,m(d%ijk_reorder)%qrs%triangular%tuple_combInd, &
           & array_tot(m(d%ijk_reorder)%qrs%triangular%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
           & m(d%ijk_reorder)%qrs%triangular%tuple,B%integral,sorted)

    else

        if (if_trans=='N') then
          row    = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%pqr%triangular%tuple_combInd
       endif

        call srt_31(save_index)%restricted_bra (nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,c(d%K)%orbital_type, &
        & c(d%L)%orbital_type,m(d%ijk_reorder)%pqr%triangular%tuple_combInd, &
        & array_tot(m(d%ijk_reorder)%pqr%triangular%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
        & m(d%ijk_reorder)%pqr%triangular%tuple,B%integral,sorted) 

     endif

      if (do_inverse) then

       if (if_trans=='N') then
          row    = c(d%bra)%triangular%tuple_combInd 
          column = c(d%ket+shift1)%tuple_combInd
       else  
          column = c(d%bra)%triangular%tuple_combInd 
          row    = c(d%ket+shift1)%tuple_combInd
       endif
       endif

         case('restricted_ket')

        call input_index(c,d,m)

        shift1 = 4

     if (save_index==0) then

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%pqr%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%pqr%tuple_combInd
       endif

     call srt_31(save_index)%restricted_ket(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,&
          & c(d%K)%orbital_type,c(d%L)%orbital_type,m(d%ijk_reorder)%pqr%tuple_combInd, &
          & array_tot(m(d%ijk_reorder)%pqr%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
          & m(d%ijk_reorder)%pqr%tuple,B%integral,sorted)

      elseif (save_index==1) then
       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%qrs%tuple_combInd
       endif

     call srt_31(save_index)%restricted_ket(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%K)%orbital_type,&
          & c(d%L)%orbital_type,m(d%ijk_reorder)%qrs%tuple_combInd, &
          & array_tot(m(d%ijk_reorder)%qrs%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
          & m(d%ijk_reorder)%qrs%tuple,B%integral,sorted)


        elseif (save_index == 3) then 


       if (if_trans=='N') then
        column = c(d%l_reorder)%orbital_type
        row = m(d%ijk_reorder)%pqr%triangular%tuple_combInd  
       else
        row = c(d%l_reorder)%orbital_type
        column = m(d%ijk_reorder)%pqr%triangular%tuple_combInd  
       endif
       
        call srt_31(save_index)%unrestricted(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%I)%orbital_type, &
           & c(d%J)%orbital_type,c(d%ket+shift1)%tuple_combInd,m(d%ijk_reorder)%pqr%tuple_combInd,&
           & array_tot(m(d%ijk_reorder)%pqr%triangular%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
           & m(d%ijk_reorder)%pqr%triangular%tuple,B%integral,sorted) 

     else

       if (if_trans=='N') then
          row    = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
          column = c(d%l_reorder)%orbital_type  
       else
          row    = c(d%l_reorder)%orbital_type  
          column = m(d%ijk_reorder)%qrs%triangular%tuple_combInd
       endif

     call srt_31(save_index)%restricted_ket(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,c(d%I)%orbital_type,&
     & c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd,m(d%ijk_reorder)%qrs%triangular%tuple_combInd,&
     & array_tot(m(d%ijk_reorder)%qrs%triangular%tuple_combInd,c(d%l_reorder)%orbital_type,nrep), &
     & m(d%ijk_reorder)%qrs%triangular%tuple,B%integral,sorted)

     endif

       if (do_inverse) then
       if (if_trans=='N') then
       row = c(d%bra+shift1)%tuple_combInd
       column = c(d%ket)%triangular%tuple_combInd     
       else 
       column = c(d%bra+shift1)%tuple_combInd
       row = c(d%ket)%triangular%tuple_combInd     
       endif
       endif

      end select

     elseif (sort_type_dummy == 13) then
      select case(arraytype)
        
         case('unrestricted')

      call input_index(c,d,m)

       shift1 = 4
       if (do_inverse) then
       if (if_trans=='N') then
        column =  c(d%ket+shift1)%tuple_combInd
        row    = c(d%bra+shift1)%tuple_combInd
       else
        row = c(d%ket+shift1)%tuple_combInd
        column = c(d%bra+shift1)%tuple_combInd
       endif
       else
       if (if_trans=='N') then
          row    = c(d%i_reorder)%orbital_type  
          column = m(d%jkl_reorder)%qrs%tuple_combInd
       else
          row    = m(d%jkl_reorder)%qrs%tuple_combInd
          column = c(d%i_reorder)%orbital_type  
       endif
       endif

      call srt_13(save_index)%unrestricted(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,&
     &  c(d%I)%orbital_type,c(d%J)%orbital_type,c(d%ket+shift1)%tuple_combInd, &
     &  array_tot(c(d%i_reorder)%orbital_type,m(d%jkl_reorder)%qrs%tuple_combInd,nrep), &
     &  m(d%jkl_reorder)%qrs%tuple,B%integral,sorted)
         
         case('restricted_braket')

       call input_index(c,d,m)

       if (do_inverse) then
       if (if_trans=='N') then
        column = c(d%ket)%triangular%tuple_combInd
        row    = c(d%bra)%triangular%tuple_combInd
       else
        row = c(d%ket)%triangular%tuple_combInd
        column = c(d%bra)%triangular%tuple_combInd
       endif
 
       else

       if (if_trans=='N') then
          row    = c(d%i_reorder)%orbital_type  
          column = m(d%jkl_reorder)%qrs%triangular%tuple_combInd
       else
          row    = m(d%jkl_reorder)%qrs%triangular%tuple_combInd
          column = c(d%i_reorder)%orbital_type  
       endif

       endif 

      call srt_13(save_index)%restricted_braket (nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,&
     &  c(d%I)%orbital_type,c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd, &
     & array_tot(c(d%i_reorder)%orbital_type,m(d%jkl_reorder)%qrs%triangular%tuple_combInd,nrep), &
     &  m(d%jkl_reorder)%qrs%triangular%tuple,B%integral,sorted)

         case('restricted_bra')

          call input_index(c,d,m)

          shift1 = 4 
      if (do_inverse) then
       if (if_trans=='N') then
          row    = c(d%bra)%triangular%tuple_combInd 
          column = c(d%ket+shift1)%tuple_combInd
       else  
          column = c(d%bra)%triangular%tuple_combInd 
          row    = c(d%ket+shift1)%tuple_combInd
       endif

       else

       if (if_trans=='N') then
          row    = c(d%i_reorder)%orbital_type  
          column = m(d%jkl_reorder)%qrs%tuple_combInd
       else
          row    = m(d%jkl_reorder)%qrs%tuple_combInd
          column = c(d%i_reorder)%orbital_type  
       endif
      endif

      call srt_13(save_index)%restricted_bra(nrep,multb,do_inverse,c(d%bra)%triangular%tuple_combInd,&
     &  c(d%I)%orbital_type,c(d%J)%orbital_type,c(d%ket+shift1)%tuple_combInd, &
     &  array_tot(c(d%i_reorder)%orbital_type,m(d%jkl_reorder)%qrs%tuple_combInd,nrep), &
     &  m(d%jkl_reorder)%qrs%tuple,B%integral,sorted)

      case('restricted_ket')

      call input_index(c,d,m)
      shift1 = 4
      if (do_inverse) then

       if (if_trans=='N') then
          row    = c(d%bra+shift1)%tuple_combInd 
          column = c(d%ket)%triangular%tuple_combInd
       else  
          column = c(d%bra+shift1)%tuple_combInd 
          row    = c(d%ket)%triangular%tuple_combInd
       endif

       else

       if (if_trans=='N') then
          row    = c(d%i_reorder)%orbital_type  
          column = m(d%jkl_reorder)%qrs%triangular%tuple_combInd
       else
          row    = m(d%jkl_reorder)%qrs%triangular%tuple_combInd
          column = c(d%i_reorder)%orbital_type  
       endif
       
      endif

      call srt_13(save_index)%restricted_ket(nrep,multb,do_inverse,c(d%bra+shift1)%tuple_combInd,&
     &  c(d%I)%orbital_type,c(d%J)%orbital_type,c(d%ket)%triangular%tuple_combInd, &
     &  array_tot(c(d%i_reorder)%orbital_type,m(d%jkl_reorder)%qrs%triangular%tuple_combInd,nrep), &
     &  m(d%jkl_reorder)%qrs%triangular%tuple,B%integral,sorted)

      end select
   endif
     if (allocated(V%presorted)) deallocate(V%presorted)
     
     if (associated(B%integral)) nullify(B%integral)
 
  end subroutine

  subroutine onebody_for_blas (array1,array2,array3,trans_1b,row,column,FreeIndx) 

!--------purpose-------------

!    arrange onebody routine for blas. get the indices for matrix multiplication as well as the end tensor.

!----------------------------


!------------------------------
   character(*),dimension(:),intent(in)  :: array1,array2,array3 ! follows the same instruction as it is in the align_array.
   character(*),intent(out)              :: trans_1b
   character(2),optional,intent(out)     :: FreeIndx(size(array1))
   integer                               :: reorder_array(size(array1))      
   integer,intent(out)                   :: row(:), column(:) 
   type(middle_indices)                  :: a(0:2)  
   integer,dimension(2)                  :: bit_reordered
   integer                               :: i
!------------------------------
#include "symm.inc"

  call align_array(array1,array2,array3,reorder_array)   

  if ((toDecimal(reorder_array)) == 21) then
  trans_1b = 'T'
  else
  trans_1b = 'N'
  endif
  
   a(0)%orbital_type => NO
   a(1)%orbital_type => NV

   bit_reordered = bitString(array1,reorder_array) 

   row = a(bit_reordered(1))%orbital_type
   column = a(bit_reordered(2))%orbital_type

   if (present(FreeIndx)) then

   FreeIndx(1:size(array1)) = array1(reorder_array(1:size(array1)))

  endif 

  end subroutine 

  subroutine swap_character(array,k,l)

    character(*),intent(inout) :: array(:)
    integer, intent(in)        :: k,l  
    character(2)               :: temp

    temp = array(k)

    array(k) = array(l)
    array(l) = temp
  end subroutine

  subroutine align_array(array1,array2,array3,reorder_array)

!-------------------------------------------------------------
    character(*),intent(in) :: array1(:), array2(:) ! array1 is the array to be reordered
    character(*),intent(in),optional :: array3(:) 
    integer, dimension(size(array1)),intent(out) :: reorder_array      ! array2 is the array from where we will extract the first indices in the reordered 
                                                                       ! array i.e, for left array it is always intermediate and for right always    
    integer                                  :: k,i                    ! left array. array3 is the array to extract rest of the indices e.g, for left array
                                                                       ! always right and for right array always intermediate.
!-------------------------------------------------------------
     k = 1
     do i = 1, size(array2)       ! loop over vorder 
        if (findloc(array1,array2(i)) > 0) then 
          reorder_array(k) = findloc(array1, array2(i))
          k=k+1
        endif
     enddo

    if (present(array3)) then 
     do i = 1, size(array3)       ! loop over vorder 
        if (findloc(array1,array3(i)) > 0) then 
          reorder_array(k) = findloc(array1,array3(i))
          k=k+1
        endif
     enddo
    endif
    
  end subroutine


! subroutine integral_fetching_out_of_core()

! bit_original = bitString (array_original, (/1,2,3,4/)) 

!  value = 10
!  if (popcnt(bit_original)==2) then
!   value = ior(iand(bit_original(1),bit_original(2)),iand(bit_original(3),bit_original(4))) 
!  endif



!  select case(integrals) 

!   case (popcnt(bit_original)==0)
!     allocate(Voooo(NV1))
!     call getOOOO(Voooo)
!   case (popcnt(bit_original)==1)
!     allocate(Vvooo(NV2))
!     call getVOOO(Vvooo)
!   case (value==1)
!     allocate(VVVOO(NV3))
!     call getVVOO(VVVOO) 
!   case (value==0)
!     allocate(VVOVO(NV4))
!     call getVOVO(VVOVO)

!  end select

! end subroutine


  subroutine integral_fetching_incore(G)
  implicit none
!-------------------------------------------------------------
  integer :: value
  integer,dimension(4):: bit_original, bit_available
  integer :: AllocateStatus 
  type(presorting),intent(out) :: G
!-------------------------------------------------------------
#include "symm.inc"
#include "complex.inc"

    bit_original = bitString (array_original,(/1,2,3,4/)) 

    value = bin2dec(bit_original)

!here you can think of keeping all the VVVO integrals in fast memory. 

   select case(value) 

    case (0)

      allocate(G%presorted(nv1*rcw))
      call getoooo(G%presorted)

    case (1)

      allocate(G%presorted(iovoot(nrep+1)*rcw))
      call getooov(G%presorted)

    case (2)

      allocate(G%presorted(nv2*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getoovo(G%presorted) 

    case (3)

      allocate(G%presorted(nv3*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getoovv(G%presorted)

   case (4)

      allocate(G%presorted(iovoot(nrep+1)*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getovoo(G%presorted)

    case (5)

!     allocate(G%presorted(nv4*rcw),stat = AllocateStatus)
!     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!     call getovov(G%presorted)

      write(*,*)'not yet implemented'

    case (6)

!     allocate(G%presorted(nv4*rcw),stat = AllocateStatus)
!     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!     call getovvo(G%presorted)


      write(*,*)'not yet implemented'

    case (7)

      allocate(G%presorted(iovvvt(nrep+1)*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getovvv_incore(G%presorted)

!      write(*,*)'ovvv type array has to be fetched in out-of-core way.'

    case (8)

      allocate(G%presorted(nv2*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getvooo(G%presorted)

    case (9)

      allocate(G%presorted(ivoov(nrep+1)*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getvoov(G%presorted)

    case (10)

      allocate(G%presorted(nv4*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getvovo(G%presorted)

    case (11)

      allocate(G%presorted(nv5*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getvovv_incore(G%presorted)

    case (12)

      allocate(G%presorted(nv3*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getvvoo(G%presorted)

    case (14)

      allocate(G%presorted(nv5*rcw),stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      call getvvvo_incore(G%presorted)

    case default
      write(*,*)'no suitable array to fetch'
   end select

  end subroutine

  subroutine input_index(a,b,k,reorder_v)

!pointer to most of the indices available in RELCCSD module. purpose is to access them in a generalized fashion. 

!-------------------------------
 integer,dimension(:),allocatable,target :: bit_reordered, bit_original
 integer,dimension(4),intent(in),optional:: reorder_v 
! logical, intent(in)                :: do_exchange
! character(2), intent(in)           :: array_exchange(4)
 type(middle_indices),intent(out)     :: a(0:30)  
 type(bit_processing),intent(out)   :: b 
 type(more_indices),intent(out),optional :: k(0:15)
 integer :: i
 integer, pointer :: bitpqr(:),bitqrs(:),bitpqr_reorder(:),bitqrs_reorder(:)
!-------------------------------
#include "symm.inc"
     a(0)%tuple => JJOO         
     a(1)%tuple => JJOV
     a(2)%tuple => JJVO
     a(3)%tuple => JJVV

     a(4)%tuple => IIOOT         
     a(5)%tuple => IIOV
     a(6)%tuple => IIVO
     a(7)%tuple => IIVV

     a(0)%tuple_combInd => MOO
     a(1)%tuple_combInd => MOV
     a(2)%tuple_combInd => MVO
     a(3)%tuple_combInd => MVV

     a(4)%tuple_combInd => NOO
     a(5)%tuple_combInd => NOV
     a(6)%tuple_combInd => NVO
     a(7)%tuple_combInd => NVV

     a(0)%triangular%tuple_combInd => NOOT
     a(3)%triangular%tuple_combInd => NVVT

!dimensions for the arrays where both bra and ket are totally symmetric.


      a(0)%full_array => JOOOO
      a(1)%full_array => JOOOV
      a(2)%full_array => JOOVO
      a(3)%full_array => JOOVV
      a(4)%full_array => JOVOO
      a(5)%full_array => JOVOV
      a(6)%full_array => JOVVO
      a(7)%full_array => JOVVV
      a(8)%full_array => JVOOO
      a(12)%full_array => JVVOO
      a(13)%full_array => JVVOV
      a(10)%full_array => JVOVO
      a(11)%full_array => JVOVV
      a(9)%full_array => JVOOV
      a(14)%full_array => JVVVOI4

      a(15)%full_array => J2OOOO     
      a(16)%full_array => J2OOOV
      a(17)%full_array => J2OOVO
      a(18)%full_array => J2OOVV
      a(19)%full_array => J2OVOO
      a(20)%full_array => J2OVOV
      a(21)%full_array => J2OVVO
      a(22)%full_array => J2OVVV
      a(23)%full_array => J2VOOO
      a(27)%full_array => J2VVOO
      a(28)%full_array => J2VVOV
      a(25)%full_array => J2VOVO
      a(26)%full_array => J2VOVV
      a(24)%full_array => J2VOOV
      a(29)%full_array => J2VVVO


!    a(4)%full_array = INT(JVVVO,4) !WARNING!!! JVVVOT is int*8 type. I cast it as int*4 from int*8

!dimensions for the arrays where both bra and ket are not necessarily totally symmetric.

!dimension for the single orbital index
     a(0)%orbital_type => NO
     a(1)%orbital_type => NV

  if (present(k)) then

     k(4)%pqr%tuple => KKVOO
     k(5)%pqr%tuple => KKVOV
     k(2)%pqr%tuple => KKOVO
     k(3)%pqr%tuple => KKOVV

     k(0)%pqr%triangular%tuple => KKOOOT
     k(1)%pqr%triangular%tuple => KKOOVT
     k(6)%pqr%triangular%tuple => KKVVOT
     k(7)%pqr%triangular%tuple => KKVVVT

     k(2)%pqr%tuple_combInd => NOVO2
     k(3)%pqr%tuple_combInd => NOVV
     k(4)%pqr%tuple_combInd => NVOO
     k(5)%pqr%tuple_combInd => NVOV

     k(0)%pqr%triangular%tuple_combInd => NOOOT
     k(1)%pqr%triangular%tuple_combInd => NOOVT
     k(6)%pqr%triangular%tuple_combInd => NVVOT
     k(7)%pqr%triangular%tuple_combInd => NVVVT

     k(1)%qrs%tuple => LLOOV
     k(2)%qrs%tuple => LLOVO
     k(6)%qrs%tuple => LLVVO
     k(5)%qrs%tuple => LLVOV
     k(0)%qrs%triangular%tuple => LLOOOT
     k(4)%qrs%triangular%tuple => LLVOOT
     k(3)%qrs%triangular%tuple => LLOVVT
     k(7)%qrs%triangular%tuple => LLVVVT

     k(2)%qrs%tuple_combInd => NOVO
     k(1)%qrs%tuple_combInd => NOOV
     k(6)%qrs%tuple_combInd => NVVO
     k(5)%qrs%tuple_combInd => NVOV2
     k(0)%qrs%triangular%tuple_combInd => NOOOT2
     k(4)%qrs%triangular%tuple_combInd => NVOOT
     k(3)%qrs%triangular%tuple_combInd => NOVVT
     k(7)%qrs%triangular%tuple_combInd => NVVVT2

     k(0)%pqr%full_array => KOOOOT
     k(1)%pqr%full_array => KOOOVT
     k(2)%pqr%full_array => KOOVOT
     k(3)%pqr%full_array => KOOVVT
     k(4)%pqr%full_array => KOVOO
     k(8)%pqr%full_array => KVOOO
     k(9)%pqr%full_array => KOVOV
     k(5)%pqr%full_array => KVOOV
     k(6)%pqr%full_array => KOVVOT
     k(11)%pqr%full_array=> KVOVV
     k(10)%pqr%full_array=> KVOVO
     k(12)%pqr%full_array=> KVVOOT
     k(13)%pqr%full_array=> KVVOVT
!    k(14)%pqr%full_array=  INT(KVVVOT,4)  !WARNING!!! KVVVOT is int*8 type. I cast it as int*4 from int*8
      
     k(0)%qrs%full_array => LOOOOT
     k(1)%qrs%full_array => LVOOOT
     k(2)%qrs%full_array => LOOVO
     k(3)%qrs%full_array => LOOVVT
     k(4)%qrs%full_array => LOVOOT
     k(5)%qrs%full_array => LOVOV
     k(6)%qrs%full_array => LOVVO 
     k(8)%qrs%full_array => LVOOOT 
     k(9)%qrs%full_array => LVOOVT 
     k(10)%qrs%full_array=> LVOVO 
     k(11)%qrs%full_array=> LVOVVT 
     k(12)%qrs%full_array=> LVVOOT
     k(14)%qrs%full_array=> LVVVO


     k(0)%totsym%full_array => IOOOOTT
!    k(2)%totsym%full_array => IOOVO
!    k(3)%totsym%full_array => IOOVVT
     k(4)%totsym%full_array => IOVOOT
     k(6)%totsym%full_array => IOVVO 
     k(7)%totsym%full_array => IOVVVT
     k(8)%totsym%full_array => IVOOOT 
     k(9)%totsym%full_array => IVOOV  
     k(10)%totsym%full_array=> IVOVO 
!    k(11)%totsym%full_array=> IVOVVT 
     k(12)%totsym%full_array=> IVVOOT
!    k(14)%totsym%full_array=> IVVVO

  endif

  if (present(reorder_v)) reorder_left = reorder_v  

  allocate (bit_reordered(size(array_original)))
  allocate (bit_original(size(array_original)))

  if (do_exchange) then
  bit_reordered = bitString (exchange,reorder_left)
  bit_original  = bitString (exchange, (/1,2,3,4/)) 
  else
  bit_reordered  = bitString (array_original,reorder_left)
  bit_original  = bitString (array_original, (/1,2,3,4/)) 
  endif

! unique indexing for pair indices.  

    b%bra_reorder = ior(bit_reordered(1),bit_reordered(2)) + iand(bit_reordered(1),bit_reordered(2)) &
    &                + ior(bit_reordered(1),iand(bit_reordered(1),bit_reordered(2))) 
    b%ket_reorder = ior(bit_reordered(3),bit_reordered(4)) + iand(bit_reordered(3),bit_reordered(4)) &
    &                + ior(bit_reordered(3),iand(bit_reordered(3),bit_reordered(4)))  

    b%braket = bin2dec(bit_reordered)


    b%bra = ior(bit_original(1),bit_original(2)) + iand(bit_original(1),bit_original(2)) &
    &        + ior(bit_original(1),iand(bit_original(1),bit_original(2)))
    b%ket = ior(bit_original(3),bit_original(4)) + iand(bit_original(3),bit_original(4)) &
    &        + ior(bit_original(3),iand(bit_original(3),bit_original(4)))

! unique indexing for combined triple indices.

    bitpqr=>bit_original(1:3)
    bitqrs=>bit_original(2:4)

    b%ijk = bin2dec(bitpqr) 
    b%jkl = bin2dec(bitqrs)

    bitpqr_reorder=>bit_reordered(1:3)  
    bitqrs_reorder=>bit_reordered(2:4)  

    b%ijk_reorder = bin2dec(bitpqr_reorder) 
    b%jkl_reorder = bin2dec(bitqrs_reorder)

! unique indexing for individual label of an array.

    b%I = bit_original(1) 
    b%J = bit_original(2)
    b%K = bit_original(3)
    b%L = bit_original(4)

    b%I_reorder = bit_reordered(1) 
    b%J_reorder = bit_reordered(2)
    b%K_reorder = bit_reordered(3)
    b%L_reorder = bit_reordered(4)
! unique indexing for full array when triangular index tuple are being used.
    b%ijkl_reorder = bin2dec(bit_reordered)  
    b%ijkl = bin2dec(bit_original)  

  deallocate (bit_reordered)
  deallocate (bit_original)
  end subroutine

  integer function findloc(array,value)

  character(*), dimension(:),intent(in) :: array
  character(2),intent(in) :: value
  integer      :: i  

     do i = 1, size(array)       ! loop over vorder 
        if (array(i)==value) then
          findloc = i
          exit
          else
          findloc = 0
        endif
     enddo
   end function

  function bitString (generic_array,reorder_index_array)

  character(*), dimension(:),intent(in) :: generic_array
  integer,      dimension(:),intent(in) :: reorder_index_array
  integer      :: i,k
  integer, dimension(size(generic_array)) ::bitString 
  character(2) :: temp

     do i = 1, size(generic_array)
        temp=generic_array(reorder_index_array(i))
        k = ichar(temp(1:1))-ichar('o')
        bitString(i)=k
     enddo

  end function   

  integer function toDecimal(index_array,trans)
   
  integer, dimension(:),intent(in),target :: index_array
  logical, optional, intent(in) :: trans
  integer :: scratch(size(index_array))
  integer,pointer :: scratch_bra(:)
  integer,pointer :: scratch_ket(:)
  logical                       :: do_trans 
  integer                       :: temp, i
  
  scratch = index_array
  if (present(trans)) then
   do_trans = trans
  else
   do_trans = .false.
  endif
 
  if (do_trans) then

! it should not be case dependent. if possible try to modify it.

  select case (sortType)
     case (22)
     scratch_bra => index_array(3:4) 
     scratch_ket => index_array(1:2) 
     case (31)
     scratch_bra => index_array(4:4)
     scratch_ket => index_array(1:3)
     case (13)
     scratch_bra => index_array(2:4)
     scratch_ket => index_array(1:1)
  end select  

     scratch = (/scratch_bra,scratch_ket/)  
  
  endif   

    toDecimal=0

    do i = size(index_array), 1 , -1  
      toDecimal = toDecimal + scratch(i)*10**(size(index_array)-i)
    enddo 

  end function

  integer function bin2dec(bitstring)

    integer, intent(in) :: bitstring(:)
    integer             :: i 

    bin2dec=0
 
    do i = 1, size(bitstring)
     bin2dec = bin2dec  + bitstring(i)*2**(size(bitstring)-i)
    enddo

  end function

  recursive subroutine permutate(E, P) 
    integer, intent(in)  :: E(:)       ! array of objects 
    integer, intent(out) :: P(:,:)     ! permutations of E 
    integer  :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
       N = size(E); Nfac = size(P,1); 
     do i=1,N                           ! cases with E(i) in front 
          if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 

         ! ... deactivate, hurt pgf90
         !forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 

         ! Brent's workaround
          do k=1, Nfac/N
            P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/)
          enddo

     end do 
  end subroutine permutate 

  function get_operator_array(character_array,integer_for_array)
  integer,intent(in) :: integer_for_array
  character(*),intent(in) :: character_array(:)  
  character(2)            :: get_operator_array(size(character_array))  
  integer :: y, temp
  integer :: x, i

  temp = integer_for_array
  i = 1

  do  
     y = temp/10**(size(character_array)-i)  
     x = mod(temp, 10**(size(character_array)-i))
     temp = x
     get_operator_array(i) = character_array(y)
     i = i+1
     if (i > size(character_array)) exit
  enddo 
  end function

  character(len=10) function seqstring(seqnum) 
    integer, intent(in) :: seqnum
    character(10)       :: temp
    write (temp,'(I0)') seqnum
    seqstring = trim(temp)
  end function

  function array_tot(row,column,nrep)

!-----------description-------------------------------------------------------
! construct the symmetry offset between an index triple bra/ket and and a single index bra/ket. 
!-----------------------------------------------------------------------------

  integer,intent(in):: nrep
  integer,intent(in):: row(nrep),column(nrep)
  integer           :: array_tot(nrep+1) 
  integer           :: irep

  array_tot(1) = 0

  do irep = 1, nrep
   array_tot(irep+1) = array_tot(irep) + row(irep)*column(irep)
  enddo

  end function

  function sorted_dimension(generic_array)

#include "symm.inc"

    character(*),intent(in) :: generic_array(:)
    integer :: sorted_dimension
    integer :: k,i,temp(4)

    k = 0 

    do i = 1, size(generic_array)
    
       temp = bitstring(generic_array,(/1,2,3,4/)) 
       if (temp(i) == 1) then
         k = k+1
       endif

    enddo   

    select case(k)

       case(0)
         
       sorted_dimension = joooo(nrep+1)  

       case(1)
         
       sorted_dimension = jooov(nrep+1)  

       case(2)

       sorted_dimension = jvovo(nrep+1)  

       case(3)

       sorted_dimension = jvvov(nrep+1) 

       case default

       write(*,*)'needs very large array, NOT an in-core case'

     end select  

  end function

  end module   
