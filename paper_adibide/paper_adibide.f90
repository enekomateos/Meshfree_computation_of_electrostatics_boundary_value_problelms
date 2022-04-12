module halton_sequence                   ! Halton sekuentzia kalkulatzeko modulua, nodoak uniformeki sakabanatzeko
 use mcf_tipos

 public :: halton

 contains
 
  function halton(b,n) result(r)         ! Funtzio honek b oinarriko halton sekuentziako lehenengo n balioak emango dizkigu bektore baten

   integer, intent(in)   :: b
   integer, intent(in)         :: n
   real(kind=dp), dimension(n) :: r
   real(kind=dp)               :: f,a
   integer                     :: i,j

   r(:)=0.0_dp                           ! Bektorea inizializatu, 0 balioarekin indize bakoitzean

   do j=1,n                              ! Indize bakoitzerako sekuentziako balioa kalkulatu
     i=j
     f=1.0_dp
     do                                  ! Do + if hau "while () do ()" egiteko da, eta ondoren algoritmoa aplikatzen da
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
 end module halton_sequence
 
 
 module funtzioak
  use mcf_tipos
  
  public  :: phi, L_ij
 
  contains
  
    function phi(i,j,epsilon)                                            ! Garapen multipolarra erabiliko dugu
     real(kind=dp), dimension(2), intent(in)    :: i, j                  ! i-k eta j-k nodoen (x,y) koordenatuak dituzte
     real(kind=dp), intent(in)                  :: epsilon
     real(kind=dp)                              :: phi, dist, r_j, r_i
     
     r_i = sqrt(i(1)**2+i(2)**2)                                         ! Nodo bakoitzaren zentroarekiko distantzia kalkulatu
     r_j = sqrt(j(1)**2+j(2)**2)
     dist = abs(r_j-r_i)                                                 ! Nodoen distantzia erlatiboa kalkulatu
     phi=sqrt(1+(epsilon*dist)**2)                                       ! Garapen mulipolarra aplikatu
    
    end function phi
    
    
    function L_ij(fi,epsilon)
     real(kind=dp), intent(in)        :: fi, epsilon
     real(kind=dp)                    :: L_ij, zatidura
     
     zatidura= (1+fi**2)/fi**3
     L_ij= zatidura*epsilon**2
     
    end function L_ij
 end module funtzioak

program paper_adibidea
 use mcf_tipos
 
 integer, parameter                 :: n=400, m=40, o= 10          ! n --> barruko nodo kopurua; m --> "boundary node" kopurua; o --> xaflako nodo kopurua
 integer                            :: i
 real(kind=dp)                      :: L, delta, r                 ! L --> xaflen luzera; delta --> xaflen y ardatzean desbiazioa zentrotik; r --> zilindroaren erradioa
 real(kind=dp), dimension(n)        :: x_bek, y_bek
 real(kind=dp), dimension(o)        :: xaf_pos_bek, xaf_neg_bek
 real(kind=dp), dimension(n+m+2*o)  :: b
 real(kind=dp), parameter           :: pi=acos(-1.0_dp)
  
 
end program paper_adibidea
