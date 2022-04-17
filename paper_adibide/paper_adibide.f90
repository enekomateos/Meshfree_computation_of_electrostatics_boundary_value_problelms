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
     real(kind=dp), dimension(:), intent(in)    :: i, j                  ! i-k eta j-k nodoen (x,y) koordenatuak dituzte
     real(kind=dp), intent(in)                  :: epsilon
     real(kind=dp)                              :: phi, dist, r_j, r_i
     
     r_i = sqrt(i(1)**2+i(2)**2)                                         ! Nodo bakoitzaren zentroarekiko distantzia kalkulatu
     r_j = sqrt(j(1)**2+j(2)**2)
     dist = abs(r_j-r_i)                                                 ! Nodoen distantzia erlatiboa kalkulatu
     phi=sqrt(1+(epsilon*dist)**2)                                       ! Garapen mulipolarra aplikatu
    
    end function phi
    
    
    function L_ij(phi,i,j,epsilon)
     real(kind=dp), dimension(:), intent(in)        :: i,j
     real(kind=dp), intent(in)                     :: epsilon
     real(kind=dp)                                 :: L_ij, zatidura, fi
     interface
      function phi(i,j,epsilon)
       use mcf_tipos
        real(kind=dp), dimension(:), intent(in) :: i,j
        real(kind=dp), intent(in)               :: epsilon
        real(kind=dp)                           :: phi
       end function
     end interface
     
     fi=phi(i,j)
     zatidura= (1+fi**2)/fi**3
     L_ij= zatidura*epsilon**2
     
    end function L_ij
 end module funtzioak

program paper_adibidea
 use mcf_tipos
 use halton_sequence
 use funtzioak
 
 integer, parameter                        :: n=400, m=40, o=10                ! n --> barruko nodo kopurua; m --> "boundary node" kopurua; o --> xaflako nodo kopurua
 integer                                   :: i, j
 real(kind=dp)                             :: L, delta, r, theta, fi, pos      ! L --> xaflen luzera; delta --> xaflen y ardatzean desbiazioa zentrotik; r --> zilindroaren erradioa
 real(kind=dp), dimension(n,2)             :: nodoak                           ! Nodo guztien (x,y) informazioa daukan bektorea
 real(kind=dp), dimension(m,2)             :: boundary_nodes                   ! Boundary node bakoitzaren (x,y) informazioa daukan bektorea
 real(kind=dp), dimension(o,2)             :: xaf_pos_nodo, xaf_neg_nodo       ! Xaflen nodo bakoitzaren informazioa daukan bektorea
 real(kind=dp), dimension(n+m+2*o,2)       :: guztiak
 real(kind=dp), dimension(n+m+2*o,n+m+2*o) :: A
 real(kind=dp), dimension(n+m+2*o)         :: b
 real(kind=dp), parameter                  :: pi=acos(-1.0_dp), epsilon=2.0_dp

 
 r=1.0_dp
 L= 0.7*r
 delta=0.1*r
 
 ! Barruko nodoak sortu
 open(unit=11, file="nodoak.dat", action="write", status="replace")
  nodoak(:,1)=halton(2,n)                                                                                 ! Barruko nodoen r balioak sortzeko
  nodoak(:,2)=halton(3,n)                                                                                 ! Barruko nodoen theta balioak sortzeko
  do i=1,n
   nodoak(i,2)=nodoak(i,2)*2*pi                                                                           ! theta-ren balioa [0,1]-->[0,2pi] zabaltzeko
   b(i)=0.0_dp                                                                                            ! Karga dentsitatea erdiko nodoetan 0 da.
   write(unit=11, fmt="(2f16.8)") sqrt(nodoak(i,1))*cos(nodoak(i,2)), sqrt(nodoak(i,1))*sin(nodoak(i,2))
   guztiak(i,1)=sqrt(nodoak(i,1))*cos(nodoak(i,2))                                                        ! A matrizea sortzeko nodo guztien koordenatuak batera beharko ditugu
   guztiak(i,2)=sqrt(nodoak(i,1))*sin(nodoak(i,2))                                                        ! A matrizea sortzeko nodo guztien koordenatuak batera beharko ditugu
  end do
 close(unit=11)
 
 ! Boundary nodes sortu
 open(unit=12, file="nodoak2.dat", action="write", status="replace")
  do i=1,m                                                                                                ! Boundary node-en theta angelua homogeneoki banatzeko [0,2*pi) tartean
   theta=2*pi*(i/real(m,dp))                                                                              ! Gogoratu, r=1 izango dela bounday node guztietarako
   boundary_nodes(i,1)=r*cos(theta)
   boundary_nodes(i,2)=r*sin(theta)             
   b(n+i)=0.0_dp                                                                                          ! Karga dentsitatea zilindroan 0 ezarriko dugu
   write(unit=12, fmt="(2f16.8)") boundary_nodes(i,1), boundary_nodes(i,2) 
   guztiak(n+i,1)=boundary_nodes(i,1)
   guztiak(n+i,2)=boundary_nodes(i,2)
  end do 
 close(unit=12)

 ! Xaflak sortu
 open(unit=13, file="nodoak3.dat", action="write", status="replace")
  xaf_pos_nodo(:,2)=delta                                                                                 ! Nodoen y koordenatua delta distantziara jarri zentrotik
  xaf_neg_nodo(:,2)=-delta
  do i=1,o                                                                                                ! Homogeneoki banatu x koordenatua
   pos=-L+2*l*(i/real(o,dp))
   xaf_pos_nodo(i,1)=pos
   xaf_neg_nodo(i,1)=pos
   b(n+m+i)=1.0_dp                                                                                        ! b bektorean hasierako potentziala idatzi
   b(n+m+o+i)=-1.0_dp
   write(unit=13, fmt="(2f16.8)") xaf_pos_nodo(i,1), xaf_pos_nodo(i,2)
   write(unit=13, fmt="(2f16.8)") xaf_neg_nodo(i,1), xaf_neg_nodo(i,2)
   guztiak(n+m+i,1)=xaf_pos_nodo(i,1)
   guztiak(n+m+i,2)=xaf_pos_nodo(i,2)
   guztiak(n+m+o+i,1)=xaf_neg_nodo(i,1)
   guztiak(n+m+o+i,2)=xaf_neg_nodo(i,2)
  end do
 close(unit=13) 
 
 ! A matrizea sortu

 do i=1,m+n+2*o
   do j=1,m+n+2*o
     if (i<n+1) then
        A(i,j)=L_ij(phi,guztiak(i,:),guztiak(j,:),20.0_dp)
     else
        A(i,j)=phi(guztiak(i,:),guztiak(j,:),20.0_dp)
     end if
   end do
 end do

 ! Sistema ebatzi behar dugu orain
 call gaussj(A,b)                                                        ! moduluak intent(inout) itxura dauka beraz gure soluzioa b matrizea izango da
 
end program paper_adibidea
