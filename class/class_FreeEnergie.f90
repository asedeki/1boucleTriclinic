module class_FreeEnergie
  use parameters
  use class_interaction
  use class_outputs
  implicit none
  type,public:: FreeEnergie
     real(kind = wp)::E
     real(kind = wp)::dE
     integer::Neq
   contains
     procedure::FE_pack
     procedure::getValues
     procedure::FE_Derivative
     procedure::FEwrite
  end type FreeEnergie
  interface FreeEnergie
     module procedure new_FreeEnergie
  end interface FreeEnergie

contains
  function new_FreeEnergie()
    type(FreeEnergie):: new_FreeEnergie
    new_FreeEnergie%E=0.0
    new_FreeEnergie%dE=0.0
    new_FreeEnergie%Neq=1
  end function new_FreeEnergie
  function FE_pack(this,ch) result(g)
    class(FreeEnergie),intent(inout):: this
    character(*),optional::ch
    real(kind=wp)::g
    if(present(ch))then
       g=this%dE
    else
       g=this%E
    end if
  end function FE_pack
  subroutine getValues(this,g)
    implicit none
    class(FreeEnergie),intent(inout)::this
    real(kind=wp)::g
!!!!!!!!!!!!!!!!!!!!
    integer::k
    this%E=g
  end subroutine getValues
  subroutine FE_Derivative(this,int,N,E0El)
    implicit none
    class(FreeEnergie),intent(out) :: this
    class(interaction),intent(in):: int
    integer,intent(in)::N
    real(kind=wp)::E0El
    integer:: i,k1,k2
    this%dE=0.0
    do k2=-N/2,N/2-1
       do k1=-N/2,N/2-1
          i=int%ind(k1,k2,k1)%i
          this%dE=this%dE+(2*int%g(i)%g2-int%g(i)%g1)
       end do
    end do
    this%dE=E0El*this%dE/real(N**2)/8.0
  end subroutine FE_Derivative
  subroutine FEwrite(this,tp2,T,out)
    class(FreeEnergie)::this
    real(kind=wp),intent(in)::tp2,T
    type(Outputs),intent(in)::out
    write(out%FreeE,'(f7.4,2x,e15.8,2x,e15.6)') &
                  tp2,T,this%E
  end subroutine FEwrite
end module class_FreeEnergie
