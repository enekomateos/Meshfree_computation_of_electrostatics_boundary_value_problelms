module halton_sequence                   ! Halton sekuentzia kalkulatzeko modulua, nodoak uniformeki sakabanatzeko
 use mcf_tipos

 public :: halton

 contains
 
  function halton(b,n) result(r)         ! Funtzio honek b oinarriko halton sekuentziako lehenengo n balioak emango dizkigu bektore baten

   integer, intent(in)         :: b
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
     real(kind=dp), dimension(:), intent(in)    :: i, j                  ! i-k eta j-k nodoen (x,y) koordenatuak dituzte
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
 use mcf_slineales
 use halton_sequence
 use funtzioak
 
 integer, parameter                    :: n=400, m=40, o=10                ! n --> barruko nodo kopurua; m --> "boundary node" kopurua; o --> xaflako nodo kopurua
 integer                               :: i, j
 real(kind=dp)                         :: L, delta, r, pos, fi             ! L --> xaflen luzera; delta --> xaflen y ardatzean desbiazioa zentrotik; r --> zilindroaren erradioa
 real(kind=dp), dimension(n,2)         :: nodoak                           ! Nodo guztien (x,y) informazioa daukan bektorea
 real(kind=dp), dimension(m,2)         :: boundary_nodes                   ! Boundary node bakoitzaren (r,theta) informazioa daukan bektorea
 real(kind=dp), dimension(o,2)         :: xaf_pos_nodo, xaf_neg_nodo       ! Xaflen nodo bakoitzaren informazioa daukan bektorea
 real(kind=dp), dimension(n+m+o,n+m+o) :: A
 real(kind=dp), dimension(n+m+2*o)     :: b
 real(kind=dp), parameter              :: pi=acos(-1.0_dp), epsilon=2.0_dp
 
 r=1.0_dp
 L= 0.7*r
 delta=0.1*r
 
 ! Barruko nodoak sortu
 nodoak(:,1)=halton(2,n)                                                 ! Barruko nodoen x balioak sortzeko
 nodoak(:,2)=halton(3,n)                                                 ! Barruko nodoen y balioak sortzeko
 do i=1,n
  nodoak(i,1)=nodoak(i,1)*2-1                                            ! Nodoen balioa [-1,1] tartera zabaldu
  nodoak(i,2)=nodoak(i,2)*2-1
  b(i)=0.0_dp                                                            ! Karga dentsitatea erdiko nodoetan 0 da
 end do
 
 ! Boundary nodes sortu
 boundary_nodes(:,1)=r                                                   ! Boundary node-en r balioa beti berdina da definizioz
 do i=1,m                                                                ! Boundary node-en theta angelua homogeneoki banatzeko [0,2*pi) tartean
  pos=2*pi*(i/real(m,dp))
  boundary_nodes(i,2)=pos
  b(i+n)=0.0_dp                                                          ! Karga dentsitatea zilindroan 0 ezarriko dugu
 end do 
 
 ! Xaflak sortu
 xaf_pos_nodo(:,2)=delta                                                 ! Nodoen y koordenatua delta distantziara jarri zentrotik
 xaf_neg_nodo(:,2)=-delta
 do i=1,o                                                                ! Homogeneoki banatu x koordenatua
  pos=-L+2*l*(i/real(o,dp))
  xaf_pos_nodo(i,1)=pos
  xaf_neg_nodo(i,1)=pos
  b(n+m+i)=1.0_dp                                                        ! b bektorean hasierako potentziala idatzi
  b(n+m+o+i)=-1.0_dp
 end do
 
 ! A matrizea sortu
 do i=1, n
  do j=1, n+m+o
   if (i<=n) then
    fi=phi(nodoak(i,:),nodoak(j,:), epsilon)
   if (i>n).and.(i<=n+m) then                                                        ! Lehenengo n lerroak L-rekin bete
    A(i,j)=L_ij(fi,epsilon)
   else                                                                  ! Hurrengo guztiak phi-rekin (boundary nodes eta xaflak)
    A(i,j)=fi
   end if
  end do
 end do
 
 ! A matrizea ebatzi
 
 call gaussj(A,b)                                                        ! moduluak intent(inout) itxura dauka beraz gure soluzioa b matrizea izango da
 
end program paper_adibidea
