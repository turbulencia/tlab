PROGRAM EXTRACT 

IMPLICIT NONE  

INTEGER,        PARAMETER    :: size_element = 4
INTEGER(KIND=4),DIMENSION(6) :: domain,subdomain 
INTEGER(KIND=4),DIMENSION(3) :: stride
INTEGER(KIND=8)              :: ivar,k,hdr_size,isize_plane_i,isize_plane_o,offset_i,offset_o
INTEGER(KIND=8)              :: nxo,nyo  
INTEGER(KIND=8)              :: imode 
REAL(KIND=size_element),DIMENSION(:,:), ALLOCATABLE :: a 
CHARACTER                    :: mode*6,var


domain   =(/1,2560,  1, 840,  1,5120/)
subdomain=(/1,2560, 40,  40,  1,5120/) 
stride   =(/1,      40,       1/) 
ALLOCATE(a(domain(2),domain(4)))

DO ivar =1,1
   SELECT CASE (ivar) 
   CASE(1)
      var='1' 
   CASE(2) 
      var='2'
   CASE(3) 
      var='3'  
   END SELECT

   DO imode=1,5 
      SELECT CASE (imode)
      CASE (1)
         mode='Euu' 
      CASE (2) 
         mode='Evv'
      CASE (3)  
         mode='Eww'
      CASE (4)  
         mode='Epp'
      CASE (5)  
         mode='E11'
      END SELECT
      OPEN(42,FILE=TRIM(ADJUSTL('pow7000.'//mode)),&
           FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD') 
      OPEN(43,FILE=TRIM(ADJUSTL('pow7000.'//mode))//'.sub',&
           FORM='UNFORMATTED',ACCESS='STREAM',STATUS='UNKNOWN')
   
      hdr_size=0
      isize_plane_i=(domain(2)-domain(1)+1)*(domain(4)-domain(3)+1)
      nxo=(subdomain(2)-subdomain(1))/stride(1) +1
      nyo=(subdomain(4)-subdomain(3))/stride(2) +1
      isize_plane_o=nxo*nyo 
   
      WRITE(*,*) nxo,nyo,isize_plane_i,isize_plane_o

      DO k=subdomain(5),subdomain(6),stride(3) 
         offset_i = hdr_size + ((k-1)*isize_plane_i              ) *size_element +1
         offset_o =            (((k-1)/stride(3))*isize_plane_o  ) *size_element +1
!         WRITE(*,*) k, offset_i, offset_o, MINVAL(a),MAXVAL(a)
         READ( 42,POS=offset_i) a  
         WRITE(43,POS=offset_o) a(&
           subdomain(1):subdomain(2):stride(1),&
           subdomain(3):subdomain(4):stride(2)) 
      ENDDO
      CLOSE(42) 
      CLOSE(43)
   ENDDO
ENDDO 

END PROGRAM EXTRACT
