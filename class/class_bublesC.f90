module class_bublesC
  use parameters
  use ode
  use class_system
  use class_indice
  implicit none
  type, public :: bubleC
     real(kind = wp),allocatable::IP(:,:,:),IC(:,:,:),IP1(:,:)     
     real(kind = wp)::vec,t_perp_ini=0.0,t_perp2_ini=0.0
     real(kind = wp)::tau_perp=0.0,Phi=0.0
     real(kind = wp)::E0=0.0,Temperature=0.0,tauSurE0,Tau
     integer::N_patche
   contains
     procedure :: create
     procedure :: destruct
     procedure::values
     procedure::values_tp20
     procedure::calcul  
     procedure::Eperp
     procedure::Eperp2
  end type bubleC
  interface bubleC
     module procedure new_bubleC
  end interface bubleC
contains
  function new_bubleC(sys)
    type(system),intent(in)::sys 
    type(bubleC)::new_bubleC
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    new_bubleC%vec=2*pi/real(sys%in%N_patche)
    new_bubleC%N_patche=sys%in%N_patche
    new_bubleC%E0=sys%in%E0
    new_bubleC.Tau=sys%in%Tau
    
    call new_bubleC%create()
  end function new_bubleC
  subroutine create(this)
    class(bubleC),intent(inout) ::this
    if(.not.allocated(this%IP))then
       allocate(this%IP(-this%N_patche/2:this%N_patche/2,-this%N_patche/2:this%N_patche/2,-this%N_patche/2:this%N_patche/2))
       allocate(this%IC(-this%N_patche/2:this%N_patche/2,-this%N_patche/2:this%N_patche/2,-this%N_patche/2:this%N_patche/2))
       allocate(this%IP1(-this%N_patche/2:this%N_patche/2,0:1))
    end if
    this%IP=0.0_wp
    this%IC=0.0_wp
    this%IP1=0.0_wp
  end subroutine create
  subroutine destruct(this)
    class(bubleC),intent(inout) ::this
    deallocate(this%IP,this%IC,this%IP1)
  end subroutine destruct
  subroutine values(this,x,Temp,tpi,tp2i,tau_perp,Phi)
    class(bubleC),intent(inout) ::this
    real(kind = wp), intent(in) :: x,Temp,tpi,tp2i,tau_perp,Phi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp) E
    integer :: kp,qp,k1,sigma,i_inf,i_sup
    character(6)::typ(2)
    
    this%t_perp_ini=tpi
    this%t_perp2_ini=tp2i
    this%tau_perp=tau_perp
    this%Phi=Phi
    print*,"tau_perp,Phi=",tau_perp,Phi
    this%Temperature=Temp

    
    typ(1)="C";typ(2)="P"
    E=this%E0
    this%E0=this%E0*exp(-x)
    this.tauSurE0=this.Tau*exp(x)/E
    
    this%IC=0.0_wp
    this%IP=0.0_wp
    this%IP1=0.0_wp
    i_inf=-this%N_patche/2
    i_sup=this%N_patche/2!-1
    !call timestamp()
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(qp) 
!$OMP DO
    do qp=i_inf,i_sup
       this%IP1(qp,0)=this%calcul(typ(2),qp,0,qp,"Eperp2")
       this%IP1(qp,1)=this%calcul(typ(2),qp,-i_inf,qp,"Eperp2")
       do kp = i_inf,i_sup
          do k1=i_inf,i_sup
             this%IC(k1,kp,qp)=this%calcul(typ(1),kp,qp,k1,"Eperp")
             this%IP(k1,kp,qp)=this%calcul(typ(2),kp,qp,k1,"Eperp")
             
          end do
       end do
    end do
!$OMP END DO 
!$OMP END PARALLEL     
    !print*,x,"   ",this%tauSurE0

    this%E0=E
    
    !this%IC(i_inf:i_sup,1:this%N_patche/2,-this%N_patche/2:this%N_patche/2)=&
    !     this%IC(i_sup:i_inf:-1,-1:-this%N_patche/2:-1,this%N_patche/2:-this%N_patche/2:-1)
    
    !this%IP(i_inf:i_sup,1:this%N_patche/2,-this%N_patche/2:this%N_patche/2)&
    !     =this%IP(i_sup:i_inf:-1,-1:-this%N_patche/2:-1,this%N_patche/2:-this%N_patche/2:-1)
    !print*,x,sum(this%IP(0,:,0)),sum(this%IP(1,:,0))
    !print*,x,32*this%IP(0,0,0),32*this%IC(0,0,0)

  end subroutine values
  subroutine values_tp20(this,x,Temp,tpi,tp2i,tau_perp,Phi)
    class(bubleC),intent(inout) ::this
    real(kind = wp), intent(in) :: x,Temp,tpi,tp2i,tau_perp,Phi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp) E
    integer :: kp,qp,k1,sigma,i_inf,i_sup
    character(6)::typ(2)
    
    this%t_perp_ini=tpi
    this%t_perp2_ini=tp2i
    
    this%tau_perp=tau_perp
    this%Phi=Phi
   
    this%Temperature=Temp
    typ(1)="C";typ(2)="P"
    E=this%E0
    this%E0=this%E0*exp(-x)
    this%tauSurE0=this%Tau/this%E0
    
    this%IC=0.0_wp
    this%IP=0.0_wp
    this%IP1=0.0_wp
    i_inf=-this%N_patche/2
    i_sup=this%N_patche/2!-1
 
    !call timestamp()
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(qp,kp,k1) 
!$OMP DO
    do qp=i_inf,i_sup-1
       this%IP1(qp,0)=this%calcul(typ(2),qp,0,qp,"Eperp2")
       this%IP1(qp,1)=this%calcul(typ(2),qp,-i_inf,qp,"Eperp2")
       do kp = i_inf,i_sup-1
          do k1=i_inf,i_sup-1
             this%IC(k1,kp,qp)=this%calcul(typ(1),kp,qp,k1,"Eperp")
             this%IP(k1,kp,qp)=this%calcul(typ(2),kp,qp,k1,"Eperp")
             
          end do
       end do
    end do
!$OMP END DO 
!$OMP END PARALLEL     
    !call timestamp()
    
    this%E0=E    
    
    !this%IC(i_inf:i_sup,1:i_sup,i_inf:i_sup)=&
    !     this%IC(i_sup:i_inf:-1,-1:i_inf:-1,i_sup:i_inf:-1)
    
    !this%IP(i_inf:i_sup,i_inf:i_sup,i_inf:0)=&
    !     this%IC(i_sup:i_inf:-1,i_sup:i_inf:-1,0:i_sup)
    !this%IP(i_inf:i_sup,i_inf:i_sup,0:i_sup)=&
    !     this%IC(i_sup:i_inf:-1,i_sup:i_inf:-1,i_inf:0)
    

  end subroutine values_tp20
  function calcul(this,typ,kp,qp,k1,eperptyp)
    class(bubleC),intent(inout) ::this
    character(*)::typ
    integer,intent(in) ::kp,qp,k1
    character(*)::eperptyp
    real(kind = wp) :: calcul
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp) :: err=1E-7,abserr=1E-7
    integer,parameter::neq=1
    real ( kind = wp ) work(100+21*neq)
    integer :: iwork(5)
    integer:: iflag
    real (kind = wp) :: k_perp, q_perp
    real(kind=wp) :: k_perp_start, k_perp_end,y(neq),ytmp(neq)
    real(kind=wp) ::k1_perp

    
    k_perp =real(kp);k1_perp =real(k1); q_perp =real(qp)
    k_perp_start = k_perp-0.5_wp ; k_perp_end = k_perp+0.5_wp
    y=0.0_wp
    iflag = 1
100 call ode_rkf (DI, neq, y ,  k_perp_start ,k_perp_end, err,abserr, iflag, work, iwork )
    if(iflag==4)goto 100
    !call DI(k_perp,ytmp,y) 
    if (abs(iflag) == 2)then 
       calcul= y(1)/real(this%N_patche)
    else
       print*,"Error in integration : Class Buble"
       call error(iflag)
       stop
    end if
  contains
    subroutine DI(k_perp,y1,yp) 
      
      real ( kind = wp ),intent(out) :: yp(neq)
      real(kind = wp), intent(in) ::y1(neq), k_perp
      
      real(kind = wp) :: A,arg1,arg2,res
      integer :: i,mu
      if(eperptyp=="Eperp")then
         A=this%Eperp(k_perp,q_perp,k1_perp,typ)
      else
         A=this%Eperp2(k_perp,q_perp,typ)
      end if
      mu =1
      yp = 0.0_wp
      arg1 = 0.5_wp*this%E0/this%Temperature
      do i=1,2
         arg2 = 0.5_wp*this%E0/this%Temperature + mu*0.5_wp*A/this%Temperature
         call calcul_tau(res,arg1,arg2,this%tauSurE0)
         yp = yp + theta(abs(this%E0+mu*A)-this%E0)*res
         mu = -1
      end do
    end subroutine DI
  end function calcul
  function Eperp(this,k_perp,q_perp,k1,typ)
    class(bubleC),intent(inout) ::this
    real(kind = wp),intent(in):: k_perp,q_perp,k1
    character(*)::typ
    real ( kind = wp ) :: Eperp
    !print*,"k1=",k1
    select case(typ)

    case('P')
      Eperp = this%Eperp2(k_perp,q_perp,'P') - this%Eperp2(k1,q_perp,'P')
       
    case('C')
      Eperp = this%Eperp2(k_perp,q_perp,'C') - this%Eperp2(k1,q_perp,'C')
   
    end select
  end function EPERP
  function Eperp2(this,k_perp,q_perp,typ)
    class(bubleC),intent(inout) ::this
    real(kind = wp),intent(in):: k_perp,q_perp
    character(*)::typ
    real ( kind = wp ) :: Eperp2
    ! print*,"tp2=",this%t_perp2_ini
    select case(typ)
    case('P')
      e('P')
       Eperp2 = 2*this%t_perp_ini*cos(k_perp*this%vec+this%Phi)& 
            + 2*this%t_perp2_ini*cos(2*k_perp*this%vec+2*this%Phi)&
            - 2*this%tau_perp*sin(2*k_perp*this%vec+2*this%Phi)&

            + 2*this%t_perp_ini*cos((q_perp+k_perp)*this%vec-this%Phi)&
            + 2*this%t_perp2_ini*cos(2*(q_perp+k_perp)*this%vec-2*this%Phi)&
            + 2*this%tau_perp*sin(2*(q_perp+k_perp)*this%vec-2*this%Phi)

    case('C')
      Eperp2 = 2*this%t_perp_ini*cos(k_perp*this%vec + this%Phi)& 
            + 2*this%t_perp2_ini*cos(2*k_perp*this%vec + 2*this%Phi)&
            - 2*this%tau_perp*sin(2*k_perp*this%vec + 2*this%Phi)&

            - 2*this%t_perp_ini*cos( (q_perp-k_perp)*this%vec - this%Phi)& 
            - 2*this%t_perp2_ini*cos(2*(q_perp-k_perp)*this%vec - 2*this%Phi)
            + 2*this%tau_perp*sin(2*(q_perp-k_perp)*this%vec - 2*this%Phi)

    end select
  end function EPERP2
end module class_bublesC


