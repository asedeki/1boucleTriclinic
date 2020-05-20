module class_indice
  use parameters
  use class_system
  implicit none
  private
  type, public ::indice
     integer :: i1, i2, i3 
   contains
     procedure::equal
     procedure::put
     procedure::geti4
     procedure::parite
     procedure::timeI
     procedure::chiral
     procedure::getindices
     procedure::print_index
     GENERIC :: operator(==) => equal
     GENERIC :: ASSIGNMENT(=)=> put
     
  end type indice
  type, public::index
     type(indice) :: k
     integer :: i
  end type index
  interface indice
     module procedure new_indice
  end interface indice
contains
  function new_indice(i1,i2,i3)
    integer, intent(in)::i1,i2,i3
    type(indice)::new_indice
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    new_indice%i1=i1
    new_indice%i2=i2
    new_indice%i3=i3
  end function new_indice
  logical function equal(r1, r2 ) !symetrie par changement de left<->right
    implicit none
    class(indice) , intent(in) :: r1
    type( indice ), intent( in ) :: r2
    if (r1%i1-r2%i1==0 .and.  r1%i2-r2%i2==0 .and.  r1%i3-r2%i3==0) then
       equal= .true.
    else
       equal=.false.
    end if
  end function equal
  subroutine put(this,r)
    implicit none
    class(indice),intent(inout):: this
    type(indice),intent(in):: r
    this%i1=r%i1
    this%i2=r%i2
    this%i3=r%i3
  end subroutine put
  integer function  geti4(ind,G)
    implicit none
    class(indice) , intent(in) :: ind
    integer, intent(in) :: G
    geti4=wrap2(ind%i1+ind%i2-ind%i3,G)
  end function geti4
  function parite( r ,G)
    implicit none
    type( indice) :: parite
    class( indice ), intent( in ) :: r
    integer,intent(in) :: G
    parite%i1 = wrap2(-1*r%i1,G)
    parite%i2 = wrap2(-1*r%i2,G)
    parite%i3 = wrap2(-1*r%i3,G)
  end function parite
  
  function timeI( r ,G) !double_permutation_circulaire
    implicit none
    class( indice ):: r
    integer, intent( in ) :: G
    type( indice) :: timeI
    
    timeI%i1=r%i3
    timeI%i2=r%geti4(G)
    timeI%i3=r%i1
  end function timeI
  
  function chiral( r ,G) !symetrie par changement de left<->right
    implicit none
    type( indice) :: chiral
    class( indice ), intent( in ) :: r
    integer, intent( in ) :: G
    chiral%i1 =r%i2
    chiral%i2 =r%i1
    chiral%i3 =r%geti4(G)
  end function chiral
  subroutine print_index(this,listIdx,listInd,n)
    implicit none
    class(indice),intent(inout):: this
    type(index),  allocatable,intent(in) :: listIdx(:,:,:)
    type(indice), allocatable,intent(in) :: listInd(:)
    integer,intent(in) :: n 
    integer :: i
    integer :: i1,i2,i3,i_inf, i_sup
    character(len=1024) :: filename,indices,indexs
    i_inf = -n/2 ; i_sup = n/2 -1
    indices="indices_n32.dat"
    indexs="indexes_n32.dat"
    write (filename, "(I2)")n
    indices(10:11)=trim(filename)
    indexs(10:11)=trim(filename)
  
    open (unit=11,file=indexs,action="write",status="unknown")
    open (unit=12,file=indices,action="write",status="unknown")
    write(unit=11,fmt='(a,I10)') 'total=',size(listIdx)
    do i1=i_inf,i_sup
       do i2=i_inf,i_sup
          do i3=i_inf,i_sup
             write(unit=11,fmt='(I6,2x,I3,2x,I3,2x,I3)') listIdx(i1,i2,i3)%i,i1,i2,i3
          end do
       end do
    end do
    write(unit=12,fmt='(a,I10)') 'total=',size(listInd)
    do i=1, size(listInd)
       write(unit=12,fmt=*) listInd(i)
    end do

  end subroutine print_index
  subroutine getindices(this,DiffIndices,MatInd,N_patche)
    implicit none
    class(indice),intent(inout):: this
    type(indice), allocatable ,intent(out)::DiffIndices(:)
    type(index), allocatable ,intent(out) :: MatInd(:,:,:)
    integer::N_patche

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(indice), allocatable :: G(:)
    type(indice), allocatable :: g_final(:)
    type(indice), allocatable ,target :: gg_temp(:)
    integer :: position
    type(indice) , pointer :: gg(:)
    integer i1, i2, i3 , ps, k, i_inf, i_sup
    i_inf = -N_patche/2 ; i_sup = N_patche/2 -1
    allocate(G(N_patche*N_patche*N_patche))
    allocate(MatInd(i_inf:i_sup,i_inf:i_sup,i_inf:i_sup))
    k=1 
    do i1= i_inf , i_sup
       do i2= i_inf , i_sup
          do i3= i_inf , i_sup  
             G(k) = indice(i1,i2,i3)
             MatInd(i1,i2,i3)%k=indice(50000,50000,50000)
             MatInd(i1,i2,i3)%i=50000
             k=k+1
          end do
       end do
    end do
    
    allocate(gg_temp(1))
    gg_temp(1) = indice(i_inf,i_inf,i_inf)
    MatInd(i_inf,i_inf,i_inf)%k = indice(i_inf,i_inf,i_inf)
    MatInd(i_inf,i_inf,i_inf)%i = 1
    gg => gg_temp
    ps = 2
    do i1 = 2, size(G)
       call append_inds(position,gg,g_final, G(i1),N_patche)
       if (allocated(g_final)) then
          MatInd(G(i1)%i1,G(i1)%i2,G(i1)%i3)%k = G(i1)
          MatInd(G(i1)%i1,G(i1)%i2,G(i1)%i3)%i = ps
          deallocate(gg_temp)
          allocate(gg_temp(size(g_final)))
          gg_temp(:) = g_final(:)
          gg => gg_temp
          deallocate(g_final)
          ps = ps +1
       else
          MatInd(G(i1)%i1,G(i1)%i2,G(i1)%i3)%k = gg(position)
          MatInd(G(i1)%i1,G(i1)%i2,G(i1)%i3)%i = position
       end if
       
    end do
    deallocate(G)
    allocate(DiffIndices(size(gg)))
    do i1 = 1, size(gg)
       DiffIndices(i1) = gg(i1)
    end do
    
  end subroutine getindices
  subroutine append_inds(position, list, autre, e,N_patche)
    type(indice), intent(inout), allocatable :: autre(:)
    type (indice), intent (in) :: e
    type (indice), pointer :: list(:)
    integer::N_patche
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (indice) :: val1(7)
    logical :: temp
    integer :: i
    integer, intent(out) :: position
    temp = .true.
    
    
    val1(1) = e%chiral(N_patche) ; val1(2) = e%parite(N_patche); val1(3) = e%timeI(N_patche)
    val1(4)=  val1(3)%parite(N_patche) 
    val1(5)=  val1(3)%chiral(N_patche) 
    val1(6) = val1(1)%parite(N_patche) 
    val1(7) = val1(6)%timeI(N_patche)
    
    do i=1,7
       if (test_existance(position, list, val1(i),N_patche)) then
          temp = .false.
          exit
       end if
    end do
    
    if ( temp) then
       allocate(  autre( size(list)+1) )
       autre(:size(list))=list(:)
       autre(size(list)+1)=e
    end if
  end subroutine append_inds

  function test_existance( position, list, r ,N_patche) 
    logical :: test_existance
    type (indice), pointer :: list(:)
    type( indice ), intent( in ) :: r
    integer, intent(out) :: position
    integer::N_patche
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer it,imod
    it=1
    test_existance = .false.
    do while (it <= size(list))
       imod=modulo(abs(list(it)%i1+r%i1),N_patche)+modulo(abs(list(it)%i2+r%i2),N_patche)&
            +modulo(abs(list(it)%i3+r%i3),N_patche)
       if( (list(it) == r).or.(imod==0)) then
          test_existance = .true.
          position = it
          exit 
       end if
       it=it+1
    end do
  end function test_existance
end module class_indice
