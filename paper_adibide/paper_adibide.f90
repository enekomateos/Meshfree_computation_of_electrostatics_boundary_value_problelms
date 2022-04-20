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
 
 integer, parameter                         :: n=400, m=40, o=10                ! n --> barruko nodo kopurua; m --> "boundary node" kopurua; o --> xaflako nodo kopurua
 integer                                    :: i, j, k, kon, dimen, total
 real(kind=dp)                              :: L, delta, r, theta, pos          ! L --> xaflen luzera; delta --> xaflen y ardatzean desbiazioa zentrotik; r --> zilindroaren erradioa
 real(kind=dp), dimension(n,2)              :: nodoak                           ! Erdiko nodoen (x,y) informazioa daukan bektorea
 real(kind=dp), dimension(:,:), allocatable :: A, b, guztiak                    ! Sistemaren matrizea eta hasierako baldintzak; Nodo guztiak (boundary+xafla ere) hemen daude
 real(kind=dp), parameter                   :: pi=acos(-1.0_dp), epsilon=2.0_dp
 real(kind=dp)                              :: u,x,y,c,d,f,g                    ! u--> Soluzioa puntu batean; x,y --> emaitza irudikatzeko; c,d,f,g --> emaitza irudikatzeko
 real(kind=dp), dimension(2)                :: bek
 
 r=1.0_dp
 L= 0.7*r
 delta=0.1*r
 dimen=0

 ! Barruko nodoak sortu
  open(unit=20, file="nodoak.dat", action="write", status="replace")
  nodoak(:,1)=halton(2,n)                                                                         ! Barruko nodoen r balioak sortzeko
  nodoak(:,2)=halton(3,n)                                                                         ! Barruko nodoen theta balioak sortzeko
  do i=1,n
    if (nodoak(i,1)>r-0.02_dp) then                                                               ! Nodoak ez egoteko boundarytik gertuegi
     cycle
    else
     nodoak(i,2)=nodoak(i,2)*2*pi                                                                 ! theta-ren balioa [0,1]-->[0,2pi] zabaltzeko
     x=sqrt(nodoak(i,1))*cos(nodoak(i,2))
     y=sqrt(nodoak(i,1))*sin(nodoak(i,2))
     if ((x<L+0.02_dp) .and. (((y<delta+0.02_dp) .and. (y>delta-0.02_dp)) .or. ((y>-delta+0.02_dp) .and. (y<-delta-0.02_dp)))) then
      cycle
     else
      dimen=dimen+1
      write(unit=20, fmt="(2f20.12)") x,y
     end if
    end if
  end do
  close(unit=20)
 ! Eman dimentsioak matrizeari eta guztiak bektoreari
 total=dimen+m+o*2
 allocate(A(total,total), b(total,1), guztiak(total,2))
 
 ! Barruko nodoak sartu eta b-ri balioak eman
 open(unit=20, file="nodoak.dat", action="read", status="replace")
 do i=1, dimen
  b(i,1)=0.0_dp                                                                                    ! Karga dentsitatea erdiko nodoetan 0 da.
  read(unit=20, fmt="(2f20.12)") guztiak(i,1), guztiak(i,2)
 end do
 close(unit=20)
 
 ! Boundary nodes sortu
  do i=1,m                                                                                         ! Boundary node-en theta angelua homogeneoki banatzeko [0,2*pi) tartean
   theta=2*pi*(i/real(m,dp))                                                                       ! Gogoratu, r=1 izango dela boundary node guztietarako
   guztiak(dimen+i,1)=r*cos(theta)
   guztiak(dimen+i,2)=r*sin(theta)             
   b(dimen+i,1)=0.0_dp                                                                                 ! Karga dentsitatea zilindroan 0 ezarriko dugu
  end do 

 ! Xaflak sortu
  do i=1,o                                                                                         ! Homogeneoki banatu x koordenatua
   pos=-L+2*l*(i/real(o,dp))
   guztiak(dimen+m+i,1)=pos
   guztiak(dimen+m+i,2)=delta
   guztiak(dimen+m+o+i,1)=pos
   guztiak(dimen+m+o+i,2)=-delta
   b(dimen+m+i,1)=1.0_dp                                                                               ! b bektorean hasierako potentziala idatzi
   b(dimen+m+o+i,1)=-1.0_dp
  end do
 close(unit=13) 
 
 ! A matrizea sortu
 do i=1,total
   do j=1,total
     if (i<dimen+1) then
        A(i,j)=L_ij(phi,guztiak(i,:),guztiak(j,:),20.0_dp)
     else
        A(i,j)=phi(guztiak(i,:),guztiak(j,:),20.0_dp)
     end if
   end do
 end do


 ! Sistema ebatzi behar dugu orain
 call gaussj(A,b)                                                                                 ! moduluak intent(inout) itxura dauka beraz gure soluzioa b matrizea izango da

 ! Ekuazio diferentziala ebatzi dugunez irudikatu dezagun emaitza  
 c=-1.0_dp
 d=1.0_dp
 f=-1.0_dp
 g=1.0_dp
 
  open(unit=11, status="replace", action="write", file="paper_datuak.dat")
  do i=1,20
     x=f+(i-1)/real(20-1)*(g-f)
     bek(1)=x
     do k=1,20
     	y=c+(k-1)/real(20-1)*(d-c)
        bek(2)=y
        if (x**2+y**2<1) then
           u=0.0_dp
           do j=1,dimen
           u=u+b(j,1)*phi(bek,guztiak(j,:),20.0_dp)
           end do
           write(unit=11, fmt="(3f20.10)") x, y, u
        else
           cycle
        end if 
     end do
  end do 

end program paper_adibidea
