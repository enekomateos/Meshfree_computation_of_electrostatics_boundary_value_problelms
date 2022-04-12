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
