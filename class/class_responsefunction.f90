module class_susc
  use parameters
  implicit none

  type,public:: susc
     real(kind = wp),allocatable :: z(:,:)
     real(kind = wp),allocatable :: chi(:)
     integer::Neq
     integer::dimz1
     integer::dimz2
   contains
     procedure::vec
     procedure::get
     procedure::eq
     procedure::Sallocate
     procedure::deallocate
     GENERIC :: ASSIGNMENT(=) => eq
  end type susc
  interface susc
     module procedure new_susc
  end interface susc

contains
  subroutine Sallocate(this)
    class(susc)::this
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(this%chi(1:this%dimz1))
    allocate(this%z(1:this%dimz1,-this%dimz2/2:this%dimz2/2))

  end subroutine Sallocate
  function new_susc(dchi,N,valini,ty)
    integer,intent(in)::dchi
    integer,intent(in)::N
    real(wp),intent(in):: valini
    character(4),optional::ty
    type(susc)::new_susc
    !!!!!!!!!!!!!!
    integer::i,j,i_inf,i_sup
    real(wp)::vec
    vec=2*pi/real(N)
    i_inf=-N/2
    i_sup=N/2

    new_susc%dimz1=dchi;new_susc%dimz2=N+1
    new_susc%Neq=dchi*(N+1)+dchi
    
    call new_susc%Sallocate()
    new_susc%chi=0.0
    if(valini==1.0)then
       if(present(ty))then
          select case(ty)
          case('1sin')
             do j=1,dchi
                new_susc%z(j,:)=sqrt(2.0_wp)*sin((/(vec*real(i),i=i_inf,i_sup)/))
             end do
          case('1cos')
             do j=1,dchi
                new_susc%z(j,:)=sqrt(2.0_wp)*cos((/(vec*real(i),i=i_inf,i_sup)/))
             end do
          case('3cos')
             do j=1,dchi
                new_susc%z(j,:)=sqrt(2.0_wp)*cos((/(3.0*vec*real(i),i=i_inf,i_sup)/))
             end do
          case('2sin')
             do j=1,dchi
                new_susc%z(j,:)=sqrt(2.0_wp)*sin((/(2.0*vec*real(i),i=i_inf,i_sup)/))
             end do
          case('2cos')
             do j=1,dchi
                new_susc%z(j,:)=sqrt(2.0_wp)*cos((/(2.0*vec*real(i),i=i_inf,i_sup)/))
             end do
          end select
       else
          new_susc%z=1.0
       end if
    else
       new_susc%z=0.0
    end if
    
  end function new_susc
  subroutine eq(this,chi)
    implicit none
    class(susc),intent(out):: this
    type(susc),intent(in):: chi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print*,"this=",allocated(this%z)
    !print*,"chi=",allocated(chi%z)
    allocate(this%chi(1:chi%dimz1))
    allocate(this%z(1:chi%dimz1,-chi%dimz2/2:chi%dimz2/2))
    
    this%z(:,:)=chi%z(:,:)
    this%chi(:)=chi%chi(:)
    this%dimz1=chi%dimz1
    this%dimz2=chi%dimz2
    this%Neq=chi%Neq
  end subroutine eq
  subroutine deallocate(this)
    class(susc),intent(inout):: this
    deallocate(this%z)
    deallocate(this%chi)
  end subroutine deallocate
  function vec(this) result(g)
    class(susc),intent(inout):: this
    real(kind=wp),allocatable::g(:)
    integer::i,j,k,N
    allocate(g(this%Neq))
    !k=size(this%chi)
    g(1:this%dimz1)=this%chi(:)
    g(this%dimz1+1:this%Neq)=reshape(this%z,(/this%Neq-this%dimz1/))
  end function vec
  subroutine get(this,g)
    implicit none
    class(susc),intent(inout)::this
    real(kind=wp)::g(:)
!!!!!!!!!!!!!!!!!!!!
    integer::k
    !k=size(this%chi)
    this%chi(:)=g(1:this%dimz1)
    this%z(:,:)=reshape(g(this%dimz1+1:this%Neq),(/this%dimz1,this%dimz2/))
  end subroutine get
end module class_susc
module class_ResponseFunctions
  use parameters
  use class_system
  use class_susc
  use class_interaction
  use class_bublesC
  use class_outputs
  use class_indice

  implicit none
!  private
  type, public ::  ResponseFunction
     integer::N_patche
     integer::Neq
     type(susc):: CSDW!Charge Site Densite wave
     type(susc):: CBDW!Charge Bond Densite wave
     type(susc):: SSDW!Spin Site Densite wave : 1-> x,2->y,3->z
     type(susc):: SBDW!Spin Bond Densite wave x
     type(susc):: SupraTrip(4)!1->px, 2->py; 3->h; 4->f
     type(susc):: SupraSing(5) !1--> s,2->dxy, 3->dx2y2; 4->g; 5->i
   contains

     procedure::eqr
     procedure::RFallocate
     procedure::RFdeallocate
     procedure::RF_pack
     procedure::getValues
     procedure::SupraDerivativeC
     procedure::RFDerivativeC
     procedure::RFwrite
     procedure::RFwriteMax
     GENERIC :: ASSIGNMENT(=) => eqr
  end type ResponseFunction

  interface ResponseFunction
     module procedure new_ResponseFunction
  end interface ResponseFunction
  
contains
  subroutine RFallocate(this,chi)
    implicit none
    class(ResponseFunction),intent(out):: this
    type(ResponseFunction),intent(in):: chi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::i
    call this%CSDW%Sallocate()
    call this%CBDW%Sallocate()
    call this%SSDW%Sallocate()
    call this%SBDW%Sallocate()
    do i=1,4
       call this%SupraTrip(i)%Sallocate()
       call this%SupraSing(i)%Sallocate()
    end do
    call this%SupraSing(5)%Sallocate()
  end subroutine RFallocate
  subroutine eqr(this,chi)
    implicit none
    class(ResponseFunction),intent(out):: this
    type(ResponseFunction),intent(in):: chi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    this%N_patche=chi%N_patche
    this%Neq=chi%Neq
    call this%RFallocate(chi)
    this%CSDW=chi%CSDW
    this%CBDW=chi%CBDW
    this%SSDW=chi%SSDW
    this%SBDW=chi%SBDW
    this%SupraTrip(:)=chi%SupraTrip(:)
    this%SupraSing(:)=chi%SupraSing(:) 
  end subroutine eqr
  function  new_ResponseFunction(N_patche,valini)
    integer,intent(in)::N_patche
    real(wp),intent(in)::valini
    type(ResponseFunction)::new_ResponseFunction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    integer::i,d
    integer::N,dim1
    N=N_patche
    dim1=2
    new_ResponseFunction%N_patche=N
    new_ResponseFunction%CSDW=susc(dim1,N,valini)
    new_ResponseFunction%CBDW=susc(dim1,N,valini)
    new_ResponseFunction%SSDW=susc(dim1,N,valini)
    new_ResponseFunction%SBDW=susc(dim1,N,valini)
    dim1=1
    !SupraTrip(4)!1->px, 2->py; 3->h; 4->f
    !SupraSing(5) !1--> s,2->dxy, 3->dx2y2; 4->g; 5->i
    new_ResponseFunction%SupraTrip(1)=susc(dim1,N,valini)
    new_ResponseFunction%SupraTrip(2)=susc(dim1,N,valini,"1sin")
    new_ResponseFunction%SupraTrip(3)=susc(dim1,N,valini,"2cos")
    new_ResponseFunction%SupraTrip(4)=susc(dim1,N,valini,"1cos")
    
    new_ResponseFunction%SupraSing(1)=susc(dim1,N,valini)
    new_ResponseFunction%SupraSing(2)=susc(dim1,N,valini,"1sin")
    new_ResponseFunction%SupraSing(3)=susc(dim1,N,valini,"1cos")
    new_ResponseFunction%SupraSing(4)=susc(dim1,N,valini,"2sin")
    new_ResponseFunction%SupraSing(5)=susc(dim1,N,valini,"3cos")
    !do i=-N/2,N/2
    !   print*,i,new_ResponseFunction%SupraSing(1)%z(1,i),new_ResponseFunction%SupraSing(2)%z(1,i)
    !end do
    !stop
    new_ResponseFunction%Neq=2*new_ResponseFunction%CSDW%Neq+2*new_ResponseFunction%SSDW%Neq
    do i=1,4
       new_ResponseFunction%Neq=new_ResponseFunction%Neq+new_ResponseFunction%SupraTrip(i)%Neq
    end do
    do i=1,5
       new_ResponseFunction%Neq=new_ResponseFunction%Neq+new_ResponseFunction%SupraSing(i)%Neq
    end do
    
  end function new_ResponseFunction
  subroutine RFdeallocate(this)
    class(ResponseFunction),intent(inout) :: this
    integer::i

    call this%CSDW%deallocate()
    call this%CBDW%deallocate()
    call this%SSDW%deallocate()
    call this%SBDW%deallocate()
    do i=1,4
       call this%SupraTrip(i)%deallocate()
    end do
    do i=1,5
       call this%SupraSing(i)%deallocate()
    end do
  end subroutine RFdeallocate
  function RF_pack(this) result(g)
    class(ResponseFunction),intent(inout):: this
    real(kind=wp),allocatable::g(:)
    integer::i,j,k,N
    N=this%Neq
    allocate(g(N))
    k=1
    g(k:this%CSDW%Neq)=this%CSDW%vec();k=k+this%CSDW%Neq
    g(k:k+this%CBDW%Neq-1)=this%CBDW%vec();k=k+this%CBDW%Neq
    
    g(k:k+this%SSDW%Neq-1)=this%SSDW%vec();k=k+this%SSDW%Neq
    g(k:k+this%SBDW%Neq-1)=this%SBDW%vec();k=k+this%SBDW%Neq
    
    do i=1,4
       g(k:k+this%SupraTrip(i)%Neq-1)=this%SupraTrip(i)%vec()
       k=k+this%SupraTrip(i)%Neq
    end do
    do i=1,5
       g(k:k+this%SupraSing(i)%Neq-1)=this%SupraSing(i)%vec()
       k=k+this%SupraSing(i)%Neq
    end do
  end function RF_pack
  subroutine getValues(this,g)
    implicit none
    class(ResponseFunction),intent(inout):: this
    real(kind=wp)::g(:)
    integer::i,k=0
    call this%CSDW%get(g(1:this%CSDW%Neq));k=this%CSDW%Neq
    call this%CBDW%get(g(k+1:k+this%CBDW%Neq));k=k+this%CBDW%Neq
    
    call this%SSDW%get(g(k+1:k+this%SSDW%Neq));k=k+this%SSDW%Neq
    call this%SBDW%get(g(k+1:k+this%SBDW%Neq));k=k+this%SBDW%Neq
    
    do i=1,4
       call this%SupraTrip(i)%get(g(k+1:k+this%SupraTrip(i)%Neq))
       k=k+this%SupraTrip(i)%Neq
    end do
    do i=1,5
       call this%SupraSing(i)%get(g(k+1:k+this%SupraSing(i)%Neq))
       k=k+this%SupraSing(i)%Neq
    end do
    
  end subroutine getValues 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Correction!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SupraDerivativeC(this,ksi,b,g)
    implicit none
    class(ResponseFunction):: this
    type(ResponseFunction),intent(in) :: ksi
    type(bubleC),intent(in):: b
    class(interaction),intent(in):: g
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    integer::i,i_inf,i_sup,qperp(2)
    
    i_inf=-ksi%N_patche/2
    i_sup=ksi%N_patche/2-1
    qperp=0
    qperp(1)=0;qperp(2)=1!ksi%N_patche/2
    
    do i=1,2
       this%CSDW%chi(i)=dot_product((ksi%CSDW%z(i,i_inf:i_sup))**2,b%IP1(i_inf:i_sup,qperp(i)))
       this%CBDW%chi(i)=dot_product((ksi%CBDW%z(i,i_inf:i_sup))**2,b%IP1(i_inf:i_sup,qperp(i)))
       this%SSDW%chi(i)=dot_product((ksi%SSDW%z(i,i_inf:i_sup))**2,b%IP1(i_inf:i_sup,qperp(i)))
       this%SBDW%chi(i)=dot_product((ksi%SBDW%z(i,i_inf:i_sup))**2,b%IP1(i_inf:i_sup,qperp(i)))
    end do
    
    do i=1,4
       this%SupraTrip(i)%chi(1)=dot_product((ksi%SupraTrip(i)%z(1,i_inf:i_sup))**2,b%IC(0,i_inf:i_sup,qperp(1)))
    end do
    do i=1,5
       this%SupraSing(i)%chi(1)=dot_product((ksi%SupraSing(i)%z(1,i_inf:i_sup))**2,b%IC(0,i_inf:i_sup,qperp(1)))
    end do

  end subroutine SupraDerivativeC
  subroutine RFDerivativeC(this,ksi,b,g)
    implicit none
    class(ResponseFunction):: this
    type(ResponseFunction),intent(in) :: ksi
    type(bubleC),intent(in):: b
    class(interaction),intent(in):: g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer:: i,j,kpp
    integer::sigma,qperp(2),i_inf,i_sup,kperp1
    real(wp)::Ipc,Ipp,Ic,kc,kp
    integer::ind1,ind2,ind3,ind4,ind5
    

    call this%SupraDerivativeC(ksi,b,g)
    

    i_inf=-ksi%N_patche/2
    i_sup=ksi%N_patche/2-1
    !kperp1=-i_inf/2

    
    kc=0.0
    kp=-i_inf

    this%CSDW%z=0.0;this%CBDW%z=0.0
    this%SSDW%z=0.0;this%SBDW%z=0.0
    do j=1,4
       this%SupraTrip(j)%z=0.0
       this%SupraSing(j)%z=0.0
    end do
    this%SupraSing(5)%z=0.0
    do kpp = i_inf,i_sup
       do i =  i_inf,i_sup
          ind1=g%ind(i,kpp,kpp)%i
          ind2=g%ind(i,kpp-kp,kpp)%i
          ind3=g%ind(i,kpp,i-kp)%i
          ind4=g%ind(i,kpp,kpp-kp)%i
          ind5=g%ind(kpp,-kpp,i)%i
          
          Ipc=b%IP(kpp,i,kc)
          Ipp=b%IP(kpp,i,kp)
          Ic=b%IC(kpp,i,kc)

         this%CSDW%z(1,kpp) =this%CSDW%z(1,kpp) -0.5*(2*g%g(ind1)%g1-g%g(ind1)%g2+&!effet g3
               2*g%g(ind1)%g3-g%g(g%ind(i,kpp,i)%i)%g3)*ksi%CSDW%z(1,i)*Ipc

         this%CSDW%z(2,kpp) =this%CSDW%z(2,kpp) -0.5*(2*g%g(ind2)%g1-g%g(ind2)%g2+&!effet g3
               2*g%g(ind4)%g3-g%g(ind3)%g3)*ksi%CSDW%z(2,i)*Ipp

         this%CBDW%z(1,kpp) =this%CBDW%z(1,kpp) -0.5*(2*g%g(ind1)%g1-g%g(ind1)%g2+&!effet g3
               -2*g%g(ind1)%g3+g%g(g%ind(i,kpp,i)%i)%g3)*ksi%CBDW%z(1,i)*Ipc            
         this%CBDW%z(2,kpp) =this%CBDW%z(2,kpp) -0.5*(2*g%g(ind2)%g1-g%g(ind2)%g2+&!effet g3
               -2*g%g(ind4)%g3+g%g(ind3)%g3)*ksi%CBDW%z(2,i)*Ipp

          
         this%SSDW%z(1,kpp) =this%SSDW%z(1,kpp) +0.5*(g%g(ind1)%g2+&!effet g3
               g%g(g%ind(i,kpp,i)%i)%g3)*ksi%SSDW%z(1,i)*Ipc
          
         this%SSDW%z(2,kpp) =this%SSDW%z(2,kpp) +0.5*(g%g(ind2)%g2+&!effet g3
               g%g(ind3)%g3)*ksi%SSDW%z(2,i)*Ipp
          
          
         this%SBDW%z(1,kpp) =this%SBDW%z(1,kpp) +0.5*(g%g(ind1)%g2+&!effet g3
               -g%g(g%ind(i,kpp,i)%i)%g3)*ksi%SBDW%z(1,i)*Ipc
          
         this%SBDW%z(2,kpp) =this%SBDW%z(2,kpp) +0.5*(g%g(ind2)%g2+&!effet g3
               -g%g(ind3)%g3)*ksi%SBDW%z(2,i)*Ipp
         
         !type(susc):: SupraTrip(4)!1->px, 2->py; 3->h; 4->f
         !type(susc):: SupraSing(5) !1--> s,2->dxy, 3->dx2y2; 4->g; 5->i
         do j=1,4
            this%SupraTrip(j)%z(1,kpp)  = this%SupraTrip(j)%z(1,kpp) + &
                 0.5*(g%g(ind5)%g1-g%g(ind5)%g2)*ksi%SupraTrip(j)%z(1,i)*Ic
            this%SupraSing(j)%z(1,kpp)  = this%SupraSing(j)%z(1,kpp) - &
                 0.5*(g%g(ind5)%g1+g%g(ind5)%g2)*ksi%SupraSing(j)%z(1,i)*Ic   
         end do
         this%SupraSing(5)%z(1,kpp)  = this%SupraSing(5)%z(1,kpp) - &
              0.5*(g%g(ind5)%g1+g%g(ind5)%g2)*ksi%SupraSing(5)%z(1,i)*Ic
      end do
   end do
 end subroutine RFDerivativeC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End Correction !!!!!!!!!!!!!!!!!!!!!
 subroutine RFwrite(this,tp2,T,out)
    class(ResponseFunction)::this
    real(kind=wp),intent(in)::tp2,T
    type(Outputs),intent(in)::out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::i

    write(out%chiCDW,'(f7.4,2x,e15.8,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)') &
         tp2,T,this%CSDW%chi(1),this%CBDW%chi(1),this%CSDW%chi(2),this%CBDW%chi(2)

    write(out%chiSDW,'(f7.4,2x,e15.8,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)') &
         tp2,T,this%SSDW%chi(1),this%SBDW%chi(1),this%SSDW%chi(2),this%SBDW%chi(2)

    write(out%chiST,'(f7.4,2x,e15.8,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)') &
         tp2,T,this%SupraTrip(1)%chi(1),this%SupraTrip(2)%chi(1),this%SupraTrip(3)%chi(1),this%SupraTrip(4)%chi(1)

    write(out%chiSS,'(f7.4,2x,e15.8,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)') &
         tp2,T,this%SupraSing(1)%chi(1),this%SupraSing(2)%chi(1),this%SupraSing(3)%chi(1),this%SupraSing(4)%chi(1),this%SupraSing(5)%chi(1)
    call this%RFwriteMax(tp2,T,out)
  end subroutine RFwrite
  subroutine RFwriteMax(this,tp2,T,out)
    class(ResponseFunction)::this
    real(kind=wp),intent(in)::tp2,T
    type(Outputs),intent(in)::out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::i
    real(kind=wp)::maxCDW,maxSDW,maxST,maxSS
    maxCDW=max(this%CSDW%chi(1),this%CBDW%chi(1),this%CSDW%chi(2),this%CBDW%chi(2))
    maxSDW=max(this%SSDW%chi(1),this%SBDW%chi(1),this%SSDW%chi(2),this%SBDW%chi(2))
    
    maxST=max(this%SupraTrip(1)%chi(1),this%SupraTrip(2)%chi(1),this%SupraTrip(3)%chi(1),this%SupraTrip(4)%chi(1))
    maxSS=max(this%SupraSing(1)%chi(1),this%SupraSing(2)%chi(1),this%SupraSing(3)%chi(1),this%SupraSing(4)%chi(1),&
         this%SupraSing(5)%chi(1))
    
    write(out%chiMax,'(f7.4,2x,e15.8,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)') &
         tp2,T,maxCDW,maxSDW,maxST,maxSS
  end subroutine RFwriteMax
end module class_ResponseFunctions
