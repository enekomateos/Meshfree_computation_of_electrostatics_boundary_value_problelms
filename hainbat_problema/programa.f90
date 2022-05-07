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
     real(kind=dp)                              :: phi, dist, r_x, r_y, r_z
     
     r_x=(i(1)-j(1))**2                                                  ! x direkzioko distantzia
     r_y=(i(2)-j(2))**2                                                  ! y direkzioko distantzia
     r_z=(i(3)-j(3))**2
     dist=r_x+r_y+r_z
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
 
 integer, parameter                        :: n=1000, m=110                 ! n --> barruko nodo kopurua; m --> "boundary node" kopurua; o --> xaflako nodo kopurua
 integer                                   :: i, j, k, l, npausu, ptukop, dimen
 real(kind=dp)                             :: r1, r2, det          ! L --> xaflen luzera; delta --> xaflen y ardatzean desbiazioa zentrotik; r1 --> kanpoko geruza esferikoaren erradioa; r2 --> barneko geruzaren erradioa
 real(kind=dp), dimension(2*n,3)           :: nodoak                           ! Erdiko nodoen (x,y) informazioa daukan bektorea
 real(kind=dp), dimension(n+4*m,4)         :: guztiak                          ! Nodo guztiak (boundary+xafla ere) hemen daude
 real(kind=dp), dimension(n+4*m,n+4*m)     :: A                                ! Sistemaren matrizea
 real(kind=dp), dimension(n+4*m)           :: b           
 integer, dimension(n+4*m)                 :: indizeak                 ! Hasierako baldintzak
 real(kind=dp), parameter                  :: pi=acos(-1.0_dp), epsilon=20.0_dp
 real(kind=dp)                             :: u,x,y,z,c,d,f,g                   
 real(kind=dp), dimension(4*m)             :: x_bek,y_bek
 real(kind=dp), dimension(3)               :: bek

 r1=1.0_dp
 r2=r1/2.0_dp

 ! Barruko nodoak sortu
  nodoak(:,1)=halton(2,2*n)                                                                         ! Barruko nodoen x balioak sortzeko
  nodoak(:,2)=halton(3,2*n)                                                                         ! Barruko nodoen y balioak sortzeko
  nodoak(:,3)=halton(7,2*n)                                                                         ! Barruko nodoen z balioak sortzeko

  dimen=0
  i=0
  do 
    i=i+1
    x=nodoak(i,1)*2-1.0_dp
    y=nodoak(i,2)*2-1.0_dp
    z=nodoak(i,3)*2-1.0_dp
    if ((x**2+y**2+z**2)<r1**2) then                                                                 
     dimen=dimen+1
     guztiak(dimen,1)=x
     guztiak(dimen,2)=y
     guztiak(dimen,3)=z
     guztiak(dimen,4)=0.0_dp
     b(dimen)=0.0_dp
    end if
    if (dimen==n) then
     exit
    end if
  end do
 

 ! Boundary nodes sortu
 x_bek=halton(2,4*m)*2-1.0_dp   
 y_bek=halton(3,4*m)*2-1.0_dp   
 i=0
 j=0
 do 
     j=j+1                       
     x=x_bek(j)
     y=y_bek(j)            
     if (x**2+y**2<=r1**2) then
        i=i+1
        z=sqrt(r1**2-x**2-y**2)
        guztiak(n+i,1)=x
        guztiak(n+i,2)=y
        guztiak(n+i,3)=z   
        guztiak(n+i,4)=1.0_dp
        b(n+i)=1.0_dp
        guztiak(n+m+i,1)=x
        guztiak(n+m+i,2)=y             
        guztiak(n+m+i,3)=-z
        guztiak(n+m+i,4)=-1.0_dp             
        b(n+m+i)=-1.0_dp
     end if  
     if (i==m) then
        exit
     end if
 end do

 !BARRUKO BOUNDAY NODE-ak SORTU
 i=0
 j=0
 do
     j=j+1
     x=x_bek(j)
     y=y_bek(j)
     if (x**2+y**2<=r2**2) then
        i=i+1
        z=sqrt(r2**2-x**2-y**2)
        guztiak(n+2*m+i,1)=x
        guztiak(n+2*m+i,2)=y
        guztiak(n+2*m+i,3)=z
        guztiak(n+2*m+i,4)=-1.0_dp
        b(n+2*m+i)=-1.0_dp
        guztiak(n+3*m+i,1)=x
        guztiak(n+3*m+i,2)=y
        guztiak(n+3*m+i,3)=-z
        guztiak(n+3*m+i,4)=1.0_dp
        b(n+3*m+i)=1.0_dp
     end if
     if (i==m) then
        exit
     end if
 end do

open(unit=12, file="nodoak.dat", status="replace", action="write")
  write(unit=12,fmt="(4f22.12)") (guztiak(i,:), i=1, n+4*m)
close(unit=12)


 ! A matrizea sortu
 do i=1,n+4*m
   do j=1,n+4*m
     if (i<n+1) then
        A(i,j)=L_ij(phi,guztiak(i,:),guztiak(j,:),epsilon)
     else
        A(i,j)=phi(guztiak(i,:),guztiak(j,:),epsilon)
     end if
   end do
 end do

! Sistema ebatzi behar dugu orain
 call lu_descomposicion(A,indizeak,det)             
 call lu_resolucion(A,indizeak,b)                   

 ! Ekuazio diferentziala ebatzi dugunez irudikatu dezagun emaitza  
 c=-1.0_dp
 d=1.0_dp
 f=-1.0_dp
 g=1.0_dp
 ptukop= 200

 open(unit=11, status="replace", action="write", file="datuak.dat")
  npausu=n+4*m
!  do i=1,ptukop
!     x=c+(i-1)/real(ptukop-1)*(d-c)
     bek(1)=0.0_dp
     do j=1,ptukop
       y=f+(j-1)/real(ptukop-1)*(g-f)
       bek(2)=y
        do k=1,ptukop
           z=f+(k-1)/real(ptukop-1)*(g-f)
           bek(3)=z      
           if ((x**2+y**2+z**2)<r1**2) then
              u=0.0_dp
              do l=1,npausu
                 u=u+b(l)*phi(bek,guztiak(l,:), epsilon)
              end do
              write(unit=11, fmt="(4f30.12)") x, y, z, u
           end if 
        end do      
     end do
!  end do 
  close(unit=11)

end program paper_adibidea
