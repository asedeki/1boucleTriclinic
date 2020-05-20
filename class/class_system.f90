module class_outputs
  implicit none
  !private
  type, public :: Outputs
     character (len=64):: file_out_chiCDW='',file_out_chiSDW=''
     character (len=64):: file_out_chiST=''
     character (len=64):: file_out_chiSS=''
     character (len=64):: file_out_chiMax=''
     character (len=64):: file_out_g=''
     character (len=64):: file_out_gmoy=''
     character (len=64):: file_out_FreeE=''
     character (len=64):: file_out_resistivite=''
     character (len=64):: file_out_seebeck=''

     integer::fg=10,fgmoy=20,chiCDW=40,chiSDW=50,chiSS=100,chiST=200,chiMax=130,FreeE=140,resistivite=401,seebeck=402
   contains
     procedure::get_fileNames
  end type Outputs
  
contains
  subroutine get_fileNames(this,file_out)
    class(Outputs), intent(inout) :: this
    character(*) , intent(in) :: file_out
    
    integer:: len_fo
    len_fo=len_trim(file_out)

    this%file_out_g(1: len_fo) =file_out;this%file_out_g(len_fo+1:len_fo+6)="_g.dat"
    this%file_out_gmoy(1: len_fo) =file_out;this%file_out_gmoy(len_fo+1:len_fo+9)="_gmoy.dat"
 
    this%file_out_chiCDW(1: len_fo) =file_out;this%file_out_chiCDW(len_fo+1:len_fo+12)="_chi_CDW.dat"
    this%file_out_chiSDW(1: len_fo) =file_out;this%file_out_chiSDW(len_fo+1:len_fo+14)="_chi_SDW.dat"
    this%file_out_chiSS(1: len_fo) =file_out;this%file_out_chiSS(len_fo+1:len_fo+11)="_chi_SS.dat"
    
    this%file_out_chiST(1: len_fo) =file_out;this%file_out_chiST(len_fo+1:len_fo+11)="_chi_ST.dat"
    
    this%file_out_chiMax(1: len_fo) =file_out;this%file_out_chiMax(len_fo+1:len_fo+12)="_chi_Max.dat"
    this%file_out_FreeE(1: len_fo) =file_out;this%file_out_FreeE(len_fo+1:len_fo+16)="_FreeEnergie.dat"
    this%file_out_resistivite(1: len_fo) =file_out;this%file_out_resistivite(len_fo+1:len_fo+16)="_resistivite.dat"
    this%file_out_seebeck(1: len_fo) =file_out;this%file_out_seebeck(len_fo+1:len_fo+12)="_seebeck.dat"
    
  end subroutine get_fileNames
end module class_outputs
module class_inputs
  use parameters
  implicit none
  !private
  type, public :: Inputs
     real ( kind = wp ) :: t_perp_ini=0,t_perp2_ini=0, tau_perp=0, Phi=0
     real ( kind = wp)  :: g1_ini=0, g2_ini=0,g3_ini=0
     real ( kind = wp)  :: g1_perp_ini=0,g2_perp_ini=0,g3_perp_ini=0
     real ( kind = wp ) :: Temperature=0,Temperature_final=0,d_Temperature(4),Temperature_Precision
     real ( kind = wp ) :: E0=0
     real ( kind = wp)  :: error=0
     real ( kind = wp)  :: Tau=0.0
     
     logical            :: find_TC=.false.
     integer            :: N_patche=32
     logical            :: stop=.false.
     logical            :: correction=.false.

     
     character (len=64):: file_wn="",file_Temp=""
   contains
     procedure::read_value
  end type Inputs
  interface Intputs
     module procedure init_intputs
  end interface Intputs
contains
  type(Inputs) function init_intputs() result(this)
    this%t_perp_ini=0;this%t_perp2_ini=0
    this%tau_perp=0;this%Phi=0
    this%g1_ini=0;this%g2_ini=0;this%g3_ini=0
    this%g1_perp_ini=0;this%g2_perp_ini=0;this%g3_perp_ini=0
    this%Temperature=0;this%Temperature_final=0;this%d_Temperature=0
    this%E0=0
    this%error=0
    this%find_TC=.false.
    this%N_patche=0
    this%stop=.false.
    this%Tau=0.0_wp
  end function init_intputs
  subroutine read_value(this,name,val,ok)
    class(Inputs),intent(inout)::this
    character(len=24),intent(inout):: name,val
    logical, intent(inout)::ok
    INTEGER :: index
    this%stop=.false.
    select case(name)
    case('correction=')
       read (val,*)this%correction
       ok=.true.  
    case('Tau=')
       read (val,*)this%Tau
       ok=.true.  
    case('g1_ini=')
       read (val,*) this%g1_ini
       ok=.true.
    case('g2_ini=')
       read (val,*) this%g2_ini
       ok=.true.
    case('g1_perp_ini=')
       read (val,*) this%g1_perp_ini
       ok=.true.
    case('g2_perp_ini=')
       read (val,*) this%g2_perp_ini
       ok=.true.
    case('g3_ini=')
       read (val,*) this%g3_ini
       ok=.true.
    case('g3_perp_ini=')
       read (val,*) this%g3_perp_ini
       ok=.true.
    case('t_perp_ini=')
       read (val,*) this%t_perp_ini
       ok=.true.
    case('tau_perp=')
       read (val,*) this%tau_perp
       ok=.true.
    case('Phi=')
       read (val,*) this%Phi
       ok=.true.
    
    case('find_TC=')
       read (val,*) this%find_TC
       ok=.true.
    case('t_perp2_ini=')
       read (val,*) this%t_perp2_ini
       ok=.true.
    case('error=')
       read (val,*) this%error
       ok=.true.
    case('Temp=')
       read (val,*) this%Temperature
       ok=.true.
    case('PrecisionTemp=')
       read (val,*) this%Temperature_Precision
       ok=.true.
    case('E0=')
       read (val,*) this%E0
       ok=.true.
    case('file_Temp=')
       read(val,*) this%file_Temp
       ok=.true.
    case('dTemp=')
       read(val,*)this%d_Temperature(1)
       ok=.true.
    case('N_patche=')
       read(val,*) this%N_patche
       ok=.true.
    end select
    
  end subroutine read_value
end module class_inputs
module class_system
  use parameters
  use class_Outputs
  use class_inputs
  implicit none
  !private
  type, public :: System
     type(Outputs) :: out
     type(Inputs) :: in
     
   contains
     procedure :: read => read_parameters
     procedure :: init_file
  end type System
  interface System
     module procedure new_System
  end interface system
contains
  function new_System(fileN,ch)
    type(System)::new_System
    character(*) , intent(in) :: fileN
    character(*) , optional:: ch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(present(ch))then
       call new_system%read(fileN,ch)
    else
       call new_system%read(fileN)
    end if
  end function new_System
  subroutine read_parameters(this,fileN,ch)
    class(System), intent(inout) :: this
    character(*) , intent(in) :: fileN
    character(*) , optional:: ch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer:: len_fo
    character(len=64):: name,val,file_out,file_name
    integer(4) ::yes
    logical::ok
   
    open (unit=11,file=fileN,action="read",status="old")
    read (unit=11,fmt=*,IOSTAT=yes) name,val
     
    do while(yes==0)
       name=trim(name)
       val=trim(val)
       
       ok=.false.
       call this%in%read_value(name,val,ok)
       if(.not. ok)then
          select case(name)
          case('file_name=')
             file_name=val
          case('file_out=')
             file_out=val
          end select
       end if
       read (unit=11,fmt=*,IOSTAT=yes) name,val
    end do
    close(11)
    if(present(ch))then
       call this%out%get_fileNames(file_out)
       open (unit=this%out%fg,file=this%out%file_out_g,action="write",status="unknown")
       call this%init_file(this%out%fg,"g")

        open (unit=this%out%fgmoy,file=this%out%file_out_gmoy,action="write",status="unknown")
        call this%init_file(this%out%fgmoy,"gmoy")
       
       
       open (unit=this%out%chiCDW,file=this%out%file_out_chiCDW,action="write",status="unknown")
       call this%init_file(this%out%chiCDW,"chiCDW")
       open (unit=this%out%chiSDW,file=this%out%file_out_chiSDW,action="write",status="unknown")
       call this%init_file(this%out%chiSDW,"chiSDW")
    
       open (unit=this%out%chiST,file=this%out%file_out_chiST,action="write",status="unknown")
       call this%init_file(this%out%chiST,"chiST")
    
       open (unit=this%out%chiSS,file=this%out%file_out_chiSS,action="write",status="unknown")
       call this%init_file(this%out%chiSS,"chiSS")
       

       open (unit=this%out%chiMax,file=this%out%file_out_chiMax,action="write",status="unknown")
       call this%init_file(this%out%chiMax,"chiMax")
    
       open (unit=this%out%FreeE,file=this%out%file_out_FreeE,action="write",status="unknown")
       call this%init_file(this%out%FreeE,"FreeE")

       open (unit=this%out%resistivite,file=this%out%file_out_resistivite,action="write",status="unknown")
       call this%init_file(this%out%resistivite,"r")

       open (unit=this%out%seebeck,file=this%out%file_out_seebeck,action="write",status="unknown")
       call this%init_file(this%out%seebeck,"see")
    end if
    
  end subroutine read_parameters
  subroutine init_file(this,file,typef)
    class(System),intent(in) :: this
    integer,intent(in):: file
    character(*) ::typef
    call timestamp_file (file)
    write(file,'(a)') "##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    if(this%in%correction)then
       write(file,'(a)') "##Le programme est: Quasi1D_1B_avec_correction: "
    else
       write(file,'(a)') "##Le programme est: Quasi1D_1B_sans_correction: "
    end if
    write(file,'(a)') "##The parameters used are:"
    write(file,'(a,I4)') "##N_patche=",this%in%N_patche
    write(file,'(a, e11.4)') "##Error=", this%in%error
    write(file,'(a, e11.4)') "##E_F=", this%in%E0
    write(file,'(a,f8.2)')"##t_perp_ini=", this%in%t_perp_ini
    write(file,'(a,f5.2)')"##t_perp2_ini=", this%in%t_perp2_ini
    write(file,'(a,f5.2)')"##tau_perp=", this%in%tau_perp
    write(file,'(a,f5.2)')"##Phi=", this%in%Phi
    write(file,'(a,f5.2)')"##Tau=", this%in%Tau
    write(file,'(a,e12.4)')"##g1_ini=", this%in%g1_ini
    write(file,'(a,e12.4)')"##g2_ini=", this%in%g2_ini
    write(file,'(a,e12.4)')"##g3_ini=", this%in%g3_ini
    write(file,'(a,e11.4)')"##g1_perp_ini=", this%in%g1_perp_ini
    write(file,'(a,e11.4)')"##g2_perp_ini=", this%in%g2_perp_ini
    write(file,'(a,e11.4)')"##g3_perp_ini=", this%in%g3_perp_ini
    write(file,'(a,f10.3)')"##Temperature_ini=",this%in%Temperature
    write(file,'(a)')"##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    
    select case(typef)
    case('o')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","S1","S2","S3","S4","C1","C2","C3","C4"
    case('r')
       write(file,fmt='(a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","resistivite","Error"
    case('see')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","CoeffSeebeck","CoeffSeebeckImpurete","Error"
    case('k')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","st_f","ss_dxy","ss_dx2y2","ss_g","ss_s","ss_i","st_px","st_py","st_h","tp","tp2","mu"
    case('chiCDW')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a,3x,a)')&!,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","CSDW(q=0)","CBDW(q=0)","CSDW(q=pi)","CBDW(q=pi)"
    case('chiSDW')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","SSDW(q=0)","SBDW(q=0)","SSDW(q=pi)","SBDW(q=pi)"
    case('chiST')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","ST_px","ST_py","ST_h","ST_f"
    case('chiSS')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","SS_s","SS_dxy","SS_dx2y2","SS_g","SS_i"
    case('chiMax')
       write(file,fmt='(a,3x,a,3x,a,3x,a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","CDW","SDW","ST","SS"
    case('FreeE')
       write(file,fmt='(a,3x,a,3x,a)')&
            "##t_perp2_ini","Temp","FreeEnergie"
    case('g')
       write(file,fmt='(a,3x,a,3x,a,4x,a,4x,a,6x,a,9x,a,9x,a,9x,a)')&
            "##tp2_ini","Temp","k1","k2","k3","g1","g2","g3"
    case('s')
       write(file,fmt='(a,6x,a,6x,a,6x,a,6x,a)')&
            "##t_perp2_ini","wn","Temp","kperp","z"
    case('t')
       write(file,fmt='(a,6x,a,6x,a,6x,a)')"##t_perp2_ini","Temp","t_perp","t_perp2"
    end select
    
  end subroutine init_file  
end module class_system
