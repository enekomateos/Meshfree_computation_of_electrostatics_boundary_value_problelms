module halton_sequence
use mcf_tipos

public :: halton
contains
function halton(b,n) result(r)

   integer, intent(in)   :: b
   integer, intent(in)         :: n
   real(kind=dp), dimension(n) :: r
   real(kind=dp)               :: f,a
   integer                     :: i,j

   r(:)=0.0_dp

   do j=1,n
     i=j
     f=1.0_dp
     do
        if (i>0) then
           f=f/real(b,dp)
           A=i-floor(i/real(b))*b
           r(j)=r(j)+f*a
           i=floor(i/real(b))
        else
           exit
        end if
      end do
   end do

end function halton
end module halton_secuence

program paper_adibidea
 use mcf_tipos


  contains
  
  

end program paper_adibidea
