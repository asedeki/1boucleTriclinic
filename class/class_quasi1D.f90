module class_quasi1D
  use parameters
  use class_system
  use class_interaction
  use class_bublesC
  use class_ResponseFunctions
  use class_FreeEnergie


  implicit none
  !private
  type, public :: quasi1D
     integer::Nequations=0
     real ( kind = wp ) :: t_perp=0
     real ( kind = wp ) :: t_perp2=0
     real ( kind = wp)  :: g1_ini=0, g2_ini=0,g3_ini=0
     real ( kind = wp)  :: g1_perp_ini=0,g2_perp_ini=0,g3_perp_ini=0
     real ( kind = wp ) :: Temperature=0
     real ( kind = wp ) :: E0=0
     real ( kind = wp ) :: l=500
     real ( kind = wp ) :: lc=0,ic=0.0
     real ( kind = wp ) :: resistivity(2)
     
     real ( kind = wp ) :: TempsRelax(3)
     real ( kind = wp ) :: seebeck,seebeckImp
     
     type(System):: sys
     type(Interaction)::interaction
     type(bubleC) ::bc
     type(ResponseFunction)::RF
     type(ResponseFunction)::deriv_RF
     type(FreeEnergie):: FE
  
   contains
     procedure::Intdeallocate
     procedure::flow
     procedure::flowT
     procedure::TheDerivativeC0
     procedure::TheDerivativeC
     procedure::pack
     procedure::unpack
     procedure::WriteOutputs

     procedure :: Resistivite
     procedure :: Res
     
     procedure :: CoeffSeebeck
     procedure :: CoeffSeebeck_E
     procedure :: CoeffTau
     
     procedure::ResistiviteWrite
     procedure::SeebecktiviteWrite
     
     
     procedure::End
     procedure::equal
     
     GENERIC :: ASSIGNMENT(=)=> equal
          
  end type quasi1D
  interface quasi1D
     module procedure new_quasi1D
  end interface quasi1D
contains
  subroutine Intdeallocate(this)
    class(quasi1D),intent(out)::this
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call this%deriv_RF%RFdeallocate()
  end subroutine Intdeallocate
  subroutine equal(this,from)
    class(quasi1D),intent(out)::this
    type(quasi1D) ,intent(in)::from
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(wp)::val
    this%Nequations=from%Nequations
    this%t_perp=from%t_perp
    this%t_perp2=from%t_perp2
    this%g1_ini=from%g1_ini;this%g2_ini=from%g2_ini;this%g3_ini=from%g3_ini
    this%g1_perp_ini=from%g1_perp_ini;this%g2_perp_ini=from%g2_perp_ini
    this%g3_perp_ini=from%g3_perp_ini
    this%Temperature=from%Temperature
    this%E0=from%E0
    this%l=from%l
    this%sys=from%sys
    this%interaction=from%interaction
    this%bc=from%bc

    this%RF=from%RF
    this%deriv_RF=from%deriv_RF
    this%FE=from%FE
  end subroutine equal
  function new_quasi1D(fileN,ch)
    type(quasi1D)::new_quasi1D
    character(*) , intent(in)   :: fileN
    character(*) ,optional::ch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(wp)::val
    if(present(ch))then
       new_quasi1D%sys=System(fileN,ch)
    else
       new_quasi1D%sys=System(fileN)
    end if

    new_quasi1D%t_perp=new_quasi1D%sys%in%t_perp_ini
    new_quasi1D%t_perp2=new_quasi1D%sys%in%t_perp2_ini
    new_quasi1D%g1_ini=new_quasi1D%sys%in%g1_ini
    new_quasi1D%g2_ini=new_quasi1D%sys%in%g2_ini
    new_quasi1D%g3_ini=new_quasi1D%sys%in%g3_ini
    new_quasi1D%g1_perp_ini=new_quasi1D%sys%in%g1_perp_ini
    new_quasi1D%g2_perp_ini=new_quasi1D%sys%in%g2_perp_ini
    new_quasi1D%g3_perp_ini=new_quasi1D%sys%in%g3_perp_ini
    new_quasi1D%Temperature=new_quasi1D%sys%in%Temperature
    new_quasi1D%E0=new_quasi1D%sys%in%E0
    
    new_quasi1D%interaction=Interaction(new_quasi1D%sys)
    new_quasi1D%bc=bubleC(new_quasi1D%sys)
        
    val=1.0
    new_quasi1D%RF=ResponseFunction(new_quasi1D%sys%in%N_patche,val)
    val=0.0
    new_quasi1D%deriv_RF=ResponseFunction(new_quasi1D%sys%in%N_patche,val)
    new_quasi1D%FE=FreeEnergie()
    new_quasi1D%Nequations=new_quasi1D%interaction%Neq+new_quasi1D%RF%Neq&
         +new_quasi1D%FE%Neq
    
  end function new_quasi1D
  logical function flow(this,l_ini,l_fin)
    class(quasi1D)::this
    real(kind=wp),intent(inout)::l_ini,l_fin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=2)::cor
    if (this%t_perp2==0)then
       cor="c0"
    else
       cor="c"
    end if
    flow=this.flowT(l_ini,l_fin,cor)
    call this%Resistivite()
    call this%CoeffSeebeck()
  end function flow
  logical function flowT(this,l_ini,l_fin,cor)
    class(quasi1D)::this
    real(kind=wp),intent(inout)::l_ini,l_fin
    character(len=*),intent(in)::cor
    
    !external::Derivative,bulCalcul
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    procedure(), pointer :: bulCalcul
    procedure(), pointer :: Derivative
    integer:: flag
    real ( kind = wp ):: abs_err
    real(kind=wp),allocatable:: y(:)
    real(kind=wp),allocatable:: yp(:)
    real ( kind = wp )::RG_error,cst=0.0
    logical::arret
    select case(cor)
    case("c0")
       Derivative=>fdc0
    case("c")
       Derivative=>fdc
    end select
    flowT=.false.
    arret=.false.
    this%deriv_RF=ResponseFunction(this%sys%in%N_patche,cst)
    
    allocate(yp(this%Nequations));yp=0.0_wp
    allocate(y(this%Nequations));y=0.0_wp
    
    call this%pack(y)
    abs_err=this%sys%in%error
    RG_error=this%sys%in%error
    flag = -1
    this.ic=0.0_wp
    this.lc=0.0_wp
    baba:do while ( flag < 0 .or. flag==4)
       call rkf45_d(Derivative, this%Nequations,y,yp,l_ini,&
            l_fin,RG_error,abs_err, flag)
       if (abs(this.lc-l_ini)<=1E-6)then
          this.ic=this.ic+1.0_wp
          if (this.ic>=10)then
             arret=.true.
             this.ic=0.0_wp
             this.lc=0.0_wp
          end if
       end if
       print*,"#x=",l_ini,this.lc,this.ic,arret
       this.lc=l_ini
       if(arret) exit baba
    end do baba
    
    if (abs(flag) /= 2)then
       if(abs(flag)==6.or.arret)then
          print*,"#Critical region probably reached for Temperature=",this%Temperature
          flowT=.true.
       else
          print*,"Problem in flow function -> class_quasi1D.f90"
          call error(flag)
          stop
       end if
    else
       if(arret)then
          print*,"#Critical region probably reached for Temperature=",this%Temperature
          flowT=.true.
       else
          !print*,"#flag=",flag
          call this%unpack(y)
       end if
    end if
    deallocate(yp)
    deallocate(y)
  contains
    subroutine fdc0( t, y1, yp1 )
      implicit none
      real ( kind = wp ) :: t
      real ( kind = wp ) :: y1(this%Nequations)
      real ( kind = wp ) :: yp1(this%Nequations)
      !print*,"In fd",t,size(y1)
      call this%unpack(y1)
      call this%TheDerivativeC0(t)
      call this%pack(yp1,"dy")
    end subroutine fdc0
    subroutine fdc( t, y1, yp1 )
      implicit none
      real ( kind = wp ) :: t
      real ( kind = wp ) :: y1(this%Nequations)
      real ( kind = wp ) :: yp1(this%Nequations)
      !print*,"In fd",t,size(y1)
      call this%unpack(y1)
      call this%TheDerivativeC(t)
      call this%pack(yp1,"dy")
    end subroutine fdc
  end function flowT
  
  subroutine TheDerivativeC0(this,x)
    class(quasi1D)::this
    real(kind = wp), intent(in):: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp)::E0El
    if(abs(this%l-x)>=1E-2*this%sys%in%error)then
       call this%bc%values_tp20(x,this%Temperature,this%t_perp,this%t_perp2)
    else
       this%l=x
    end if
    !this%b%IC=0.0_wp
    call this%interaction%IntDerivativeC(this%bc)
    call this%deriv_RF%RFDerivativeC(this%RF,this%bc,this%interaction)
    E0El=exp(-x)*(1.0-exp(-x))
    call this%FE%FE_Derivative(this%interaction,this%sys%in%N_patche,E0El)
  end subroutine TheDerivativeC0

  subroutine TheDerivativeC(this,x)
    class(quasi1D)::this
    real(kind = wp), intent(in):: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = wp)::E0El
    if(abs(this%l-x)>=1E-2*this%sys%in%error)then
       call this%bc%values(x,this%Temperature,this%t_perp,this%t_perp2)
    else
       this%l=x
    end if
    !this%b%IC=0.0_wp
    call this%interaction%IntDerivativeC(this%bc)
    call this%deriv_RF%RFDerivativeC(this%RF,this%bc,this%interaction)
    E0El=exp(-x)*(1.0-exp(-x))
    call this%FE%FE_Derivative(this%interaction,this%sys%in%N_patche,E0El)
  end subroutine TheDerivativeC

  subroutine pack(this,y,ch)
    class(quasi1D)::this
    real ( kind = wp ):: y(:)
    character(*),optional::ch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if(present(ch))then
       y(1:this%interaction%Neq)=this%interaction%Int_pack(ch)
       y(this%interaction%Neq+1:this%Nequations-this%FE%Neq)=this%deriv_RF%RF_pack()
       y(this%Nequations)=this%FE%FE_pack(ch)
    else
       y(1:this%interaction%Neq)=this%interaction%Int_pack()
       y(this%interaction%Neq+1:this%Nequations-this%FE%Neq)=this%RF%RF_pack()
       y(this%Nequations)=this%FE%FE_pack()
    end if
  end subroutine pack
  subroutine unpack(this,y)
    class(quasi1D)::this
    real ( kind = wp ) :: y(:)
    integer::k
    real ( kind = wp ),allocatable::yy(:)
    
    k=this%interaction%Neq
    allocate(yy(this%interaction%Neq))
    yy=y(1:this%interaction%Neq)
    call this%interaction%Int_unpack(yy) 
    deallocate(yy)
    allocate(yy(this%RF%Neq))
    yy(:)=y(this%interaction%Neq+1:this%Nequations-this%FE%Neq)
    call this%RF%getValues(yy)
    deallocate(yy)
    call this%FE%getValues(y(this%Nequations))
  end subroutine unpack
  subroutine Resistivite(this)
    class(quasi1D):: this
!!!!!!!!!!!!!!!!!!!!!Pour l matrice !!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=wp), allocatable:: L(:,:),L2(:,:),L3(:,:)
    integer:: k1,k2,k3,k4,i,Ni,Nf,ik1,ik2
    real(kind=wp)::Vi,S1,S2,g3_1,g3_2,g3_3,constante
    integer:: j1,j2,N_patche

    N_patche=this%sys%in%N_patche
    
    !Vi=(0.001-(0.5-this%sys%in%Tau)*0.02_wp/30.0_wp)/real(N_patche)
    vi=0.00D+00!+0.0001D+00*this%sys%in%Tau/real(N_patche)
    !Vi=0.001+this%sys%in%Tau*0.02_wp/30.0_wp/real(N_patche)
    Ni=-N_patche/2
    Nf=N_patche/2-1
    !allocate(L(Ni:Nf,Ni:Nf))
    allocate(L(N_patche,N_patche),L2(N_patche,N_patche),L3(N_patche,N_patche))
    constante= 1.05*this%Temperature*Pi**2/this%E0/2.0/real(N_patche)**2
    do k1=Ni,Nf
       ik1=k1-Ni+1
       do k2=Ni,Nf
          ik2=k2-Ni+1
          L(ik1,ik2)=0.0_wp
          do k3=Ni,Nf
             do k4=Ni,Nf
                i=wrap2(k4+k3-k1,N_patche)
                S1=sigma(k1,i,k3,k4,this%Temperature,N_patche,this%t_perp,this%t_perp2)
                S2=sigma(k1,k2,k3,k4,this%Temperature,N_patche,this%t_perp,this%t_perp2)
                j1=this%interaction%ind(k1,i,k3)%i
                j2=this%interaction%ind(k1,i,k4)%i
                g3_1=abs(this%interaction%g(j1)%g3-this%interaction%g(j2)%g3)**2
                j1=this%interaction%ind(k1,k2,k3)%i
                j2=this%interaction%ind(k1,k2,k4)%i
                g3_2=abs(this%interaction%g(j1)%g3-this%interaction%g(j2)%g3)**2
                j1=this%interaction%ind(k1,k3,k2)%i
                j2=this%interaction%ind(k1,k3,k4)%i
                g3_3=abs(this%interaction%g(j1)%g3-this%interaction%g(j2)%g3)**2
                
                L(ik1,ik2)=L(ik1,ik2)+(g3_1*S1*delta(k1,k2)+g3_2*S2*delta(k1+k2,k3+k4)&
                     -2*g3_3*S2*delta(k1+k3,k2+k4))*constante                   
             end do
          end do
       end do
    end do
    this%Resistivity=0.0_wp
    call this%Res(L,N_patche)

  end subroutine Resistivite
  subroutine Res(this,L,N)
    class(quasi1D):: this
!!!!!!!!!!!!!!!!!!!!!Pour l matrice !!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=wp), dimension(N,N),intent(in):: L
    integer ( kind = 4 ), intent(in):: N
    real(kind=wp), dimension(N,N)::LL
    real(kind=wp), dimension(N,1)::BB
    INTEGER:: NMAX,LWORK
    real(kind=wp), allocatable   :: WORK(:)
    
    INTEGER(kind=4), allocatable ::  IPIV(:)
    real ( kind = wp ), allocatable :: b(:,:)
    
    integer ( kind = 4 ) info,NRHS,LDA,LDB,NN
    LL(:,:)=L(:,:)
    NMAX=N
    LWORK=64*NMAX
    allocate(WORK(LWORK),IPIV(NMAX),b(N,1))
        
    LDB=N
    LDA=N
    NN=N
    NRHS=1
    
    BB(:,1)= 1.0D+00
    b(:,1) = 1.0D+00
    CALL DSYSV("Upper",NN,NRHS,L,LDA,IPIV,b,LDB,WORK,LWORK,INFO)
    this%Resistivity(1)=real(N)/sum(b(:,1))
    this%Resistivity(2)=sum(matmul(LL,b)-BB)
    deallocate(WORK,IPIV,b)
    
  end subroutine Res




!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine CoeffSeebeck(this)
    class(quasi1D):: this
    integer:: j1,j2,N_patche
    real(kind=wp)::E,y1,y2,y3,y4
    
    E=-0.001
    call this%CoeffSeebeck_E(E)
    y1=log(this%TempsRelax(1))
    y3=log(this%TempsRelax(2))
    E=0.001
    call this%CoeffSeebeck_E(E)
    y2=log(this%TempsRelax(1))
    y4=log(this%TempsRelax(2))
    this%seebeck=-500.0_wp*(y2-y1)*Pi**3*this%Temperature/3.0
    this%seebeckImp=-500.0_wp*(y4-y3)*Pi**3*this%Temperature/3.0
  end subroutine CoeffSeebeck
  subroutine CoeffSeebeck_E(this,E)
    class(quasi1D):: this
    real(kind=wp), intent(in):: E
!!!!!!!!!!!!!!!!!!!!!Pour l matrice !!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=wp), allocatable:: L(:,:),Limp(:,:),LTot(:,:)
    integer:: k1,k2,k3,k4,i,Ni,Nf,ik1,ik2
    real(kind=wp)::Vi,S1,S2,g3_1,g3_2,g3_3,constante
    integer:: j1,j2,N_patche

    N_patche=this%sys%in%N_patche
    
    Ni=-N_patche/2
    Nf=N_patche/2-1
    allocate(L(N_patche,N_patche),Limp(N_patche,N_patche),LTot(N_patche,N_patche))
    constante= 1.05*this%Temperature*Pi**2/this%E0/2.0/real(N_patche)**2
    !constante= this%Temperature/real(N_patche)**2
    do k1=Ni,Nf
       ik1=k1-Ni+1
       do k2=Ni,Nf
          ik2=k2-Ni+1
          L(ik1,ik2)=0.0_wp
          Limp(ik1,ik2)=(1.0_wp-delta(k1,k2))*1E-3/cosh(0.5_wp*E/this%Temperature)/real(N_patche)
          do k3=Ni,Nf
             do k4=Ni,Nf
                i=wrap2(k4+k3-k1,N_patche)
                
                S1=sigmaPTE(k1,i,k3,k4,E,this%Temperature,N_patche,this%t_perp,this%t_perp2)
                S2=sigmaPTE(k1,k2,k3,k4,E,this%Temperature,N_patche,this%t_perp,this%t_perp2)
                j1=this%interaction%ind(k1,i,k3)%i
                j2=this%interaction%ind(k1,i,k4)%i
                
                g3_1=abs(this%interaction%g(j1)%g3-this%interaction%g(j2)%g3)**2
                
                j1=this%interaction%ind(k1,k2,k3)%i
                j2=this%interaction%ind(k1,k2,k4)%i
                g3_2=abs(this%interaction%g(j1)%g3-this%interaction%g(j2)%g3)**2
                
                j1=this%interaction%ind(k1,k3,k2)%i
                j2=this%interaction%ind(k1,k3,k4)%i
                g3_3=abs(this%interaction%g(j1)%g3-this%interaction%g(j2)%g3)**2
                
                L(ik1,ik2)=L(ik1,ik2)+(g3_1*S1*delta(k1,k2)&
                     +g3_2*S2*delta(k1+k2,k3+k4)&
                     -2*g3_3*S2*delta(k1+k3,k2+k4))*constante                   
             end do
          end do
          
       end do
    end do

    LTot(:,:)=L(:,:)+Limp(:,:)
    call this%CoeffTau(L,N_patche,1)
    call this%CoeffTau(LTot,N_patche,2)
  end subroutine CoeffSeebeck_E

  subroutine CoeffTau(this,L,N,iT)

    class(quasi1D):: this
!!!!!!!!!!!!!!!!!!!!!Pour l matrice !!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=wp), dimension(N,N),intent(in):: L
    integer ( kind = 4 ), intent(in):: N,iT
    real(kind=wp), dimension(N,N)::LL
    real(kind=wp), dimension(N,1)::BB
    INTEGER:: NMAX,LWORK
    real(kind=wp), allocatable   :: WORK(:)
    
    INTEGER(kind=4), allocatable ::  IPIV(:)
    real ( kind = wp ), allocatable :: b(:,:)
    
    integer ( kind = 4 ) info,NRHS,LDA,LDB,NN
    LL(:,:)=L(:,:)
    NMAX=N
    LWORK=64*NMAX
    allocate(WORK(LWORK),IPIV(NMAX),b(N,1))
        
    LDB=N
    LDA=N
    NN=N
    NRHS=1
    BB(:,1)= 1.0D+00
    b(:,1) = 1.0D+00
    CALL DSYSV("Upper",NN,NRHS,L,LDA,IPIV,b,LDB,WORK,LWORK,INFO)
    
    this%TempsRelax(iT)=sum(b(:,1))!/real(N)
    this%TempsRelax(3)=max(abs(sum(matmul(LL,b)-BB)),this%TempsRelax(3))
    deallocate(WORK,IPIV,b)
  end subroutine CoeffTau

  subroutine End(this)
    class(quasi1D)::this
    call timestamp_file(this%sys%out%fg)
    call timestamp_file(this%sys%out%chiMax)
    call timestamp_file(this%sys%out%fgmoy)
    call timestamp_file(this%sys%out%chiCDW)
    call timestamp_file(this%sys%out%chiSDW)
    call timestamp_file(this%sys%out%chiSS)
    call timestamp_file(this%sys%out%chiST)
    call timestamp_file(this%sys%out%FreeE)
    call timestamp_file(this%sys%out%resistivite)
    call timestamp_file(this%sys%out%seebeck)

  end subroutine End
  subroutine WriteOutputs(this)
    class(quasi1D)::this
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call this%interaction%InteractionWrite(this%t_perp2,this%Temperature,this%sys%out)
    call this%RF%RFwrite(this%t_perp2,this%Temperature,this%sys%out)
    call this%FE%FEwrite(this%t_perp2,this%Temperature,this%sys%out)
    call this%ResistiviteWrite()
    call this%SeebecktiviteWrite()
  end subroutine WriteOutputs
  subroutine ResistiviteWrite(this)
    class(quasi1D)::this
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(this%sys%out%resistivite,'(f8.2,2x,e15.8,2x,e11.4,2x,e11.4)')&
         this%t_perp2,this%Temperature,this%Resistivity(1),this%Resistivity(2)
  end subroutine ResistiviteWrite
  subroutine SeebecktiviteWrite(this)
    class(quasi1D)::this
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(this%sys%out%seebeck,'(f8.2,2x,e15.8,2x,e11.4,2x,e11.4,2x,e11.4)')&
         this%t_perp2,this%Temperature,this%seebeck,this%seebeckImp,this%TempsRelax(3)
  end subroutine SeebecktiviteWrite

end module class_quasi1D
