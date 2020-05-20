module class_gologie
  use parameters
  use class_system
  use class_indice
  implicit none
  private
  type,public :: gologie
     real(kind=wp):: g1,dg1
     real(kind=wp):: g2,dg2
     real(kind=wp):: g3,dg3
     type(indice) :: k
   contains
     procedure::SetGIni
     procedure::eq
     GENERIC :: ASSIGNMENT(=) => eq
  end type gologie
contains
  subroutine SetGIni(this,sys)
    implicit none
    class(gologie),intent(inout)::this
    type(system)::sys 
    integer::i_inf,i_sup,k1,k2,k3
    real(kind=wp)::vec
    i_inf=-sys%in%N_patche/2
    vec = 2*pi/real(sys%in%N_patche)
    i_sup=-i_inf-1
    this%g1 = sys%in%g1_ini+&
         2*sys%in%g1_perp_ini*cos((this%k%i3-this%k%i2)*vec)
    this%g2 = sys%in%g2_ini+&
         2*sys%in%g2_perp_ini*cos((this%k%i3-this%k%i1)*vec)
    this%g3 = sys%in%g3_ini+&
         2*sys%in%g3_perp_ini*cos((this%k%i3-this%k%i2)*vec)
    this%dg1=0.0;this%dg2=0.0;this%dg3=0.0
    !print*,sys%in%g1_ini,sum(this%g1)
  end subroutine SetGIni
  subroutine eq(g,gg)
    implicit none
    class(gologie),intent(inout):: g
    type(gologie),intent(in):: gg
    g%g1=gg%g1
    g%g2=gg%g2
    g%g3=gg%g3
    g%k=gg%k
  end subroutine eq
end module class_gologie
module class_interaction
  use parameters
  use class_system
  use class_indice
  use class_gologie
  use class_bublesC
  use class_outputs

  implicit none
  type interaction
     type(gologie),allocatable:: g(:)
     type(index), allocatable :: ind(:,:,:)
     real(kind=wp):: rho
     integer::N_patche
     integer::Neq
     
   contains
    procedure::initialization
    procedure::Int_pack
    procedure::Int_unpack
    procedure::IntDerivativeC
    procedure :: assign
    procedure :: InteractionWrite
    generic :: assignment(=) => assign
    
    !generic :: operator(+) => add
 end type interaction
 interface interaction
    module procedure new_interaction
 end interface interaction
contains
  subroutine assign(to, from)
    class(interaction), intent(out) :: to
    type(interaction), intent(in) :: from
    if(.not.allocated(to%g))then
       allocate(to%g(size(from%g)))
    end if
    if(.not.allocated(to%ind))then
       allocate(to%ind(-from%N_patche:from%N_patche,-from%N_patche:from%N_patche&
            ,-from%N_patche:from%N_patche))
    end if
    to%g(:)%g1=from%g(:)%g1
    to%g(:)%g2=from%g(:)%g2
    to%g(:)%g3=from%g(:)%g3
    to%g(:)%k=from%g(:)%k
    to%ind(:,:,:)=from%ind(:,:,:)
    to%N_patche=from%N_patche
    to%Neq=from%Neq
  end subroutine assign
  function new_interaction(sys)
    type(system),intent(in)::sys 
    type(interaction)::new_interaction
    !!!!!!!!!!!!!!
    type(index), allocatable :: tc(:,:,:)
    character(len=1024) :: filename,indices,indexs
    integer:: N,i,j,k,i_inf, i_sup 
    character (len=32):: ch
    
    new_interaction%N_patche=sys%in%N_patche
    
    indices="indices_n32.dat"
    indexs="indexes_n32.dat"
    
    write (filename, "(I2)")sys%in%N_patche
    indices(10:11)=trim(filename)
    indexs(10:11)=trim(filename)

    open (unit=11,file=indices,action="read",status="old")
    read (unit=11,fmt=*) ch,N

    new_interaction%Neq=3*N
    allocate(new_interaction%g(N))
    do i=1,N
       read (unit=11,fmt=*) new_interaction%g(i)%k
    end do
    close(11)
    i_inf = -sys%in%N_patche/2 ; i_sup = sys%in%N_patche/2 -1
    allocate(tc(i_inf:i_sup,i_inf:i_sup,i_inf:i_sup))
    open (unit=101,file=indexs,action="read",status="old")
    read (unit=101,fmt=*) ch,N
    do i=i_inf,i_sup
       do j=i_inf,i_sup
          do k=i_inf,i_sup
             read (unit=101,fmt=*) tc(i,j,k)%i,tc(i,j,k)%k
             !print*,tc(i,j,k)%k
          end do
       end do
    end do
    close(101)
    i_inf = -sys%in%N_patche; i_sup = sys%in%N_patche 
    
    allocate(new_interaction%ind(i_inf:i_sup,i_inf:i_sup,i_inf:i_sup))
    do k=i_inf,i_sup
       do j=i_inf,i_sup
          do i=i_inf,i_sup
             new_interaction%ind(i,j,k)%i=tc(wrap2(i,sys%in%N_patche),wrap2(j,sys%in%N_patche),wrap2(k,sys%in%N_patche))%i
          end do
       end do
    end do
    !print*,"dim3=",bound(new_interaction%ind, 3)
    deallocate(tc)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call new_interaction%initialization(sys)
  end function new_interaction
  subroutine initialization(this,sys)
    implicit none
    class(interaction),intent(inout)::this
    type(system),intent(in)::sys 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::sizeI,i,n
    sizeI=this%Neq/3
    do i=1,sizeI
       call this%g(i)%SetGIni(sys)
    end do
  end subroutine initialization
  
  function Int_pack(this,ch) result(y)
    implicit none
    class(interaction) :: this
    character(*),optional::ch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=wp),allocatable::y(:)
    integer:: sizeg1,k

    allocate(y(this%Neq))
    
    k=1
    sizeg1=this%Neq/3
    if(present(ch))then
       y(k:sizeg1) = this%g(:)%dg1; k=k+sizeg1
       y(k:k+sizeg1-1) = this%g(:)%dg2;k=k+sizeg1
       y(k:k+sizeg1-1) = this%g(:)%dg3;k=k+sizeg1
    else
       y(k:sizeg1) = this%g(:)%g1; k=k+sizeg1
       y(k:k+sizeg1-1) = this%g(:)%g2;k=k+sizeg1
       y(k:k+sizeg1-1) = this%g(:)%g3;k=k+sizeg1
    end if
  end function Int_pack
  subroutine Int_unpack(this,y)
    implicit none
    class(interaction) :: this
    real ( kind = wp ),intent(in):: y(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer:: sizeg1,k
    
    sizeg1=this%Neq/3
    k=1
    this%g(:)%g1=y(1:sizeg1) ; k=k+sizeg1
    this%g(:)%g2=y(k:k+sizeg1-1) ;k=k+sizeg1
    this%g(:)%g3=y(k:k+sizeg1-1) ;k=k+sizeg1
  end subroutine Int_unpack
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Correction !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine IntDerivativeC(this,bub)
    class(interaction):: this
    type(bubleC),intent(in):: bub
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    real(kind=wp)::IP,IC,IPP,IP2
    integer:: N,i1 ,kp,qp,qc,i4,qpp,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12
    type(indice):: ba
    N=this%Neq/3
    this%g(:)%dg1=0.0;this%g(:)%dg2=0.0;this%g(:)%dg3=0.0
    do i1=1, N
       ba = this%g(i1)%k
       qp=wrap2(ba%i3-ba%i2,bub%N_patche);qc=wrap2(ba%i1+ba%i2,bub%N_patche)
       qpp=wrap2(ba%i1-ba%i3,bub%N_patche)
       i4 = ba%geti4(bub%N_patche)
       do kp = -bub%N_patche/2,bub%N_patche/2-1
          
          k1=this%ind(ba%i1,ba%i2,kp)%i!tc(ba%i1,ba%i2,kp)%i
          k2=this%ind(kp,qc-kp,ba%i3)%i!tc(kp,fourI(qc-kp,N_patche),ba%i3)%i
          k3=this%ind(ba%i1,kp-qp,kp)%i!tc(ba%i1,fourI(kp-qp,N_patche),kp)%i
          k4=this%ind(kp,ba%i2,ba%i3)%i!tc(kp,ba%i2,ba%i3)%i
          k5=this%ind(ba%i1,kp,kp+qp)%i!tc(ba%i1,kp,i_1)%i
          k6=this%ind(ba%i2,kp+qp,ba%i3)%i!tc(ba%i2,i_1,ba%i3)%i
          k7=this%ind(ba%i2,kp+qp,kp)%i!tc(ba%i2,i_1,kp)%i
          k8=this%ind(kp+qp,ba%i2,ba%i3)%i!tc(i_1,ba%i2,ba%i3)%i

          k9=this%ind(ba%i2,kp+qpp,kp)%i!tc(ba%i2,i_2,kp)%i
          k10=this%ind(ba%i2,kp+qpp,i4)%i!tc(ba%i2,i_2,i4)%i
          k11=this%ind(ba%i1,kp,kp+qpp)%i!tc(ba%i1,kp,i_2)%i
          k12=this%ind(ba%i1,kp,ba%i3)%i!tc(ba%i1,kp,ba%i3)%i
          IP= bub%IP(ba%i2,kp,qp)
          IP2= bub%IP(ba%i2,kp,-qp)
          IC= bub%IC(ba%i1,kp,qc)
          IPP=bub%IP(ba%i3,kp,qpp)
          
          this%g(i1)%dg1 = this%g(i1)%dg1-0.5_wp*((this%g(k1)%g2*this%g(k2)%g1+&
               this%g(k1)%g1*this%g(k2)%g2)*IC-(this%g(k3)%g2*this%g(k4)%g1 +&
               this%g(k3)%g1*this%g(k4)%g2-2*this%g(k3)%g1*this%g(k4)%g1)*IP2)+&!effet g3
               +0.5_wp*(this%g(k5)%g3*this%g(k6)%g3 +this%g(k7)%g3*this%g(this%ind(ba%i1,kp,i4)%i)%g3-&
               2_wp*this%g(k5)%g3*this%g(k8)%g3)*IP

          this%g(i1)%dg2 = this%g(i1)%dg2+0.5_wp*(-(this%g(k1)%g2*this%g(k2)%g2+&
               this%g(k1)%g1*this%g(k2)%g1)*IC+this%g(k3)%g2*this%g(k4)%g2*IP2)+&!effet g3
               0.5_wp*this%g(this%ind(ba%i1,kp,i4)%i)%g3*this%g(k6)%g3*IP

          this%g(i1)%dg3 = this%g(i1)%dg3+0.5_wp*(this%g(k5)%g3*this%g(k7)%g2+&
               this%g(this%ind(ba%i1,kp,i4)%i)%g3*this%g(k7)%g1+this%g(k5)%g2*this%g(k7)%g3+&
               this%g(k5)%g1*this%g(k6)%g3-2*this%g(k7)%g1*this%g(k5)%g3-&
               2*this%g(k7)%g3*this%g(k5)%g1)*IP+0.5_wp*(this%g(k12)%g3*this%g(k9)%g2+&
               this%g(k10)%g3*this%g(k11)%g2)*IPP
       end do
    end do
  end subroutine IntDerivativeC
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! End Correction !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine InteractionWrite(this,tp2,T,out)
    class(Interaction)::this
    real(kind=wp),intent(in)::tp2,T
    type(Outputs),intent(in)::out
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::i1,i,j
    real(wp)::g1=0.0,g2=0.0,g3=0.0
    do i1=1,this%Neq/3
       write(out%fg,'(f7.4,2x,e15.8,2x,I3,2x,I3,2x,I3,2x,e13.4,2x,e13.4,2x,e13.4,2x,e13.4)') &
            tp2,T,this%g(i1)%k%i1,this%g(i1)%k%i2,this%g(i1)%k%i3,this%g(i1)%g1,this%g(i1)%g2,this%g(i1)%g3
    end do
    write(out%fg,*);write(out%fg,*)
    do i=-this%N_patche/2,this%N_patche/2
       do j=-this%N_patche/2,this%N_patche/2
          g1=g1+this%g(this%ind(i,-i,j)%i)%g1
          g2=g2+this%g(this%ind(i,-i,j)%i)%g2
          g3=g3+this%g(this%ind(i,-i,j)%i)%g3
       end do
    end do
    g1=g1/real((this%N_patche+1)**2)
    g2=g2/real((this%N_patche+1)**2)
    g3=g3/real((this%N_patche+1)**2)
    write(out%fgmoy,'(f7.4,2x,e15.8,2x,e13.4,2x,e13.4,2x,e13.4,2x,e13.4)') tp2,T,g1,g2,g3
    
  end subroutine InteractionWrite
end module class_interaction
