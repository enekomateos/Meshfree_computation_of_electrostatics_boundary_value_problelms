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
     
     r_i=(i(1)-j(1))**2                                                  ! x direkzioko distantzia
     r_j=(i(2)-j(2))**2                                                  ! y direkzioko distantzia
     dist=r_i+r_j
     phi=sqrt(1+dist*epsilon**2)                                         ! Garapen mulipolarra aplikatu
    
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
     
     fi=phi(i,j,epsilon)
     zatidura= (1+fi**2)/fi**3
     L_ij= zatidura*epsilon**2
     
    end function L_ij
 end module funtzioak

program paper_adibidea
 use mcf_tipos
 use halton_sequence
 use funtzioak
 use mcf_slineales
 
 integer, parameter                        :: n=2000, m=400, o=400                ! n --> barruko nodo kopurua; m --> "boundary node" kopurua; o --> xaflako nodo kopurua
 integer                                   :: i, j, k, npausu, kon, ptukop, dimen
 real(kind=dp)                             :: rb, r, theta, pos, det, dx, dy      ! r --> zilindroaren erradioa; rb --> barruko zilindroaren erradioa, d --> shifts
 real(kind=dp), dimension(n*3,2)           :: nodoak                              ! Erdiko nodoen (x,y) informazioa daukan bektorea
 real(kind=dp), dimension(n+m+o,3)         :: guztiak                             ! Nodo guztiak (boundary+zilindroa ere) hemen daude
 real(kind=dp), dimension(n+m+o,n+m+o)     :: A                                   ! Sistemaren matrizea
 real(kind=dp), dimension(n+m+o)           :: b
 integer, dimension(n+m+o)                 :: indizeak                            ! Hasierako baldintzak
 real(kind=dp), parameter                  :: pi=acos(-1.0_dp), epsilon=2.0_dp
 real(kind=dp)                             :: u,x,y,c,d,f,g                       ! u--> Soluzioa puntu batean; x,y --> emaitza irudikatzeko; c,d,f,g --> emaitza irudikatzeko
 real(kind=dp), dimension(2)               :: bek
 
 r=1.0_dp
 rb=r*0.3_dp
 
 dx=0.2_dp
 dy=0.5_dp

 ! Barruko nodoak sortu
  nodoak(:,1)=halton(2,n*3)                                                                         ! Barruko nodoen r balioak sortzeko
  nodoak(:,2)=halton(3,n*3)                                                                         ! Barruko nodoen theta balioak sortzeko

  nodoak(:,1)=nodoak(:,1)*2.0_dp-1.0_dp
  nodoak(:,2)=nodoak(:,2)*2.0_dp-1.0_dp

  dimen=0
  i=0
  do 
    i=i+1
    if (((nodoak(i,1)-dx)**2+(nodoak(i,2)-dy)**2<=rb**2).or.(r**2<=(nodoak(i,1))**2+(nodoak(i,2))**2)) then
     cycle
    else
     dimen=dimen+1
     guztiak(dimen,1)=nodoak(i,1)
     guztiak(dimen,2)=nodoak(i,2)
     guztiak(dimen,3)=1.0_dp
     b(dimen)=0.0_dp
    end if
    if (dimen==n) then
     exit
    end if
  end do
 
 ! Boundary nodes sortu
  do i=1,m                                                                                                ! Boundary node-en theta angelua homogeneoki banatzeko [0,2*pi) tartean
   theta=2*pi*(i/real(m,dp))                                                                              ! Gogoratu, r=1 izango dela boundary node guztietarako
   guztiak(n+i,1)=r*cos(theta)
   guztiak(n+i,2)=r*sin(theta)             
   guztiak(n+i,3)=2.0_dp             
   b(n+i)=0.0_dp                                                                                          ! Karga dentsitatea zilindroan 0 ezarriko dugu
  end do 

 ! Barruko zilindroa sortu
  do i=1,o                                                                                                ! Homogeneoki banatu x koordenatua
   theta=2*pi*(i/real(o,dp))
   guztiak(n+m+i,1)=rb*cos(theta)+dx
   guztiak(n+m+i,2)=rb*sin(theta)+dy
   b(n+m+i)=1.0_dp                                                                                        ! b bektorean hasierako potentziala idatzi
  end do
 close(unit=13) 

open(unit=12, file="nodoak_shifted.dat", status="replace", action="write")
   do i=1,n+m+o
      write(unit=12,fmt="(3f16.8)") guztiak(i,1), guztiak(i,2), guztiak(i,3)
   end do
close(unit=12)
 
 ! A matrizea sortu
 do i=1,m+n+o
   do j=1,m+n+o
     if (i<n+1) then
        A(i,j)=L_ij(phi,guztiak(i,:),guztiak(j,:),20.0_dp)
     else
        A(i,j)=phi(guztiak(i,:),guztiak(j,:),20.0_dp)
     end if
   end do
 end do


 ! Sistema ebatzi behar dugu orain
 call lu_descomposicion(A,indizeak,det)                                                                                          ! moduluak intent(inout) itxura dauka beraz gure soluzioa b matrizea izango da
 call lu_resolucion(A,indizeak,b)                                                                                                ! moduluak intent(inout) itxura dauka beraz gure soluzioa b matrizea izango da

 ! Ekuazio diferentziala ebatzi dugunez irudikatu dezagun emaitza  
 c=-1.0_dp
 d=1.0_dp
 
 ptukop= 100
 open(unit=11, status="replace", action="write", file="paper_datuak_shifted.dat")
  npausu=n+m+o
  do i=1,ptukop
     x=c+(i-1)/real(ptukop-1)*(d-c)
     bek(1)=x
     do k=1,ptukop
     	y=c+(k-1)/real(ptukop-1)*(d-c)
      if ((x-dx)**2+(y-dy)**2<=rb**2) then
       cycle
      end if
        bek(2)=y
        if ((x)**2+(y)**2<1) then
           u=0.0_dp
           do j=1,npausu
           u=u+b(j)*phi(bek,guztiak(j,:),20.0_dp)
           end do
           write(unit=11, fmt="(3f20.10)") x, y, u
        else
           cycle
        end if 
     end do
  end do 

end program paper_adibidea
