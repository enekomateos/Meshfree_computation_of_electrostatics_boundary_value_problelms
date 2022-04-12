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
 use halton_sequence
 
 integer, parameter                 :: n=400, m=40, o= 10          ! n --> barruko nodo kopurua; m --> "boundary node" kopurua; o --> xaflako nodo kopurua
 integer                            :: i
 real(kind=dp)                      :: L, delta, r, pos            ! L --> xaflen luzera; delta --> xaflen y ardatzean desbiazioa zentrotik; r --> zilindroaren erradioa
 real(kind=dp), dimension(n,2)      :: nodoak                      ! Nodo guztien (x,y) informazioa daukan bektorea
 real(kind=dp), dimension(n,2)      :: boundary_nodes              ! Boundary node bakoitzaren (r,theta) informazioa daukan bektorea
 real(kind=dp), dimension(o,2)      :: xaf_pos_nodo, xaf_neg_nodo  ! Xaflen nodo bakoitzaren informazioa daukan bektorea
 real(kind=dp), dimension(n+m+2*o)  :: b
 real(kind=dp), parameter           :: pi=acos(-1.0_dp)
 
 r=1.0_dp
 L= 0.7*r
 delta=0.1*r
 
 ! Barruko nodoak sortu
 nodoak(:,1)=halton(2,n)                                           ! Barruko nodoen x balioak sortzeko
 nodoak(:,2)=halton(3,n)                                           ! Barruko nodoen y balioak sortzeko
 do i=1,n
  b(i)=0.0_dp                                                      ! Karga dentsitatea erdiko nodoetan 0 da
 end do
 
 ! Boundary nodes sortu
 boundary_nodes(:,1)=r                                             ! Boundary node-en r balioa beti berdina da definizioz
 do i=1,m                                                          ! Boundary node-en theta angelua homogeneoki banatzeko [0,2*pi) tartean
  pos=2*pi*(i-1/real(m,dp))
  boundary_nodes(i,2)=pos
  b(i+n)=0.0_dp                                                    ! Karga dentsitatea zilindroan 0 ezarriko dugu
 end do 
 
 ! Xaflak sortu
 xaf_pos_nodo(:,2)=delta
 xaf_neg_nodo(:,2)=-delta
 do i=1,o
 pos=-L+2*l*(i-1/real(o-1,dp))
 xaf_pos_nodo(:,1)=pos
 xaf_neg_nodo(:,1)=pos
 b(n+m+i)=1.0_dp
 b(n+m+o+i)=-1.0_dp
 end do
 
end program paper_adibidea
