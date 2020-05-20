
program test
  use parameters
  use class_system
  use class_interaction
  use class_quasi1D
  implicit none

  character(20),parameter:: file="data.tpl"
  character(2)::ecrire
  type(quasi1D):: q
  real(kind=wp):: li=2,lf=50,Temp,x,dT,TempC
  integer::i,j
  open (unit=1,file="Tc.dat",access="append",status="unknown")
  q=quasi1D(file,ecrire)

  if(q%sys%in%find_TC)then
     call flow_Temp
  else
     call flow_one_Temp
  end if
  !call f1D
  call q%End
  open (unit=1232,file="Fini",action="write",status="unknown")
contains
  subroutine flow_Temp
    implicit none
    real(kind=wp)::erreur=0.01D+00
    integer::tj
    dT=q%sys%in%d_Temperature(1)
    erreur=q%sys%in%Temperature_Precision
    write(*,*)"#dT=",q%sys%in%d_Temperature
    Temp=q%Temperature
    do i=1,100
       li=0.0;lf=100.0
       print*,"#errer=",erreur,",  dT",dT,"  dT<= err",(dT<=erreur)
       if(.not.q%flow(li,lf))then
          call q%WriteOutputs()
          !if (q%Temperature<=10.0)then
          !   exit
          !end if
       else
          print*,"#T=",q%Temperature,",  T inf",(q%Temperature<=1E-5)
          if( dT <= erreur .or. q%Temperature<=1E-5)then
             write(1,'(f9.6,2x,e15.8)')q%t_perp2,q%Temperature
             exit
          else
             if(abs(Temp-1.0D+00)<= 1E-6)then
                erreur=0.1D+00
             end if
             TempC=Temp
             Temp=Temp+dT
             dT=dT/10.0
          end if
       end if
    
       q=quasi1D(file)
       q%Temperature=Temp
       !if(abs(dT-q%Temperature)<=1E-6)then
       !   dT=dT/10.0
       !end if
       
       q%Temperature=q%Temperature-dT
       !if(q%Temperature<=1E-5)then
       !   write(1,'(f9.6,2x,a)')q%t_perp2,"1E-10"
       !   exit
       !end if
       Temp=q%Temperature
    end do
  end subroutine flow_Temp

  subroutine flow_one_Temp
    print*,"One Temp"
    li=0.0;lf=100.0
    if(.not.q%flow(li,lf))then
       call q%WriteOutputs()
    end if
  end subroutine flow_one_Temp
  subroutine f1D
    integer::i
    
    li=0.0;lf=1
    do i=1,10
       q=quasi1D(file)
       if(.not.q%flow(li,lf))then
          print*,sum(q%interaction%g(:)%g2)!
          call q%WriteOutputs()
       end if
       li=0
       lf=lf+4
    end do
  end subroutine f1D
  
end program test
    
