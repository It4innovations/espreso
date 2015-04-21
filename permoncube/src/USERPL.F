c #define inAnsys
#ifndef inAnsys
c - in the case, functions are not called in Ansys      
c
      FUNCTION vdot(a,b,n)
c      dot product of vec. a, b     sum(a(i)*b(i))
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: a, b
      DOUBLE PRECISION tmp_r
      DOUBLE PRECISION vdot
      tmp_r = 0.0d0
      DO i = 1,n
       tmp_r=tmp_r+a(i)*b(i) 
      ENDDO
      vdot=tmp_r
      END FUNCTION 
c      
      SUBROUTINE tbuser(mat,numinp,proptb)
      IMPLICIT NONE
      INTEGER n, i
      INTEGER mat,numinp
      INTEGER promM
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(numinp) :: proptb
      DOUBLE PRECISION,  DIMENSION(6) :: ci 
      DOUBLE PRECISION,  DIMENSION(6) :: gamai 
      DOUBLE PRECISION  conk,ni0,conb1,ni,mrat,ni02,conb2,ni2,conb,Rnek 

      conk=  450.0d0  
      ni0 =  1.0d0
      conb1= 10.d0
      ni  =  1.d0
      mrat=  0.0d0
      ni02=  0.0d0
      conb2= 1.0d0
      ni2 =  0.0d0 
      conb=  10.0d0
      Rnek=  0.0d0
      ci  = (/441790.,351110.,172890.,60601.,47733.,9000./)
      gamai = (/12506.,4879.,2326.,948.,258.,50./)

      proptb(1) = conk
      proptb(2) = ni0 
      proptb(3) = conb1
      ni = proptb(4)
      do i=1,6
        proptb(2*i+3) = ci(i)
        proptb(2*i+4) = gamai(i)
      end do
      proptb(17) = mrat 
      proptb(18) = ni02 
      proptb(19) = conb2 
      proptb(20) = ni2 
      proptb(21) = conb 
      proptb(22) = Rnek 


      END SUBROUTINE
c
      SUBROUTINE vmove(v1,v2,n)
c      copy of v1 into v2
      IMPLICIT NONE
      INTEGER n, i, j
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: v1, v2 
      DO i = 1,n
        v2(i)=v1(i) 
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE vzero(a,n)
c      vector entries zeroing
      IMPLICIT NONE
      INTEGER n,i,j
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: a
      DO i = 1,n
        a(i)=0.0d0 
      ENDDO
      END SUBROUTINE
c     
      SUBROUTINE maxv(a,x,y,n,m)
c      mat vec mult  y=A*x
      IMPLICIT NONE
      INTEGER n,m,i,j
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n,m) :: a
      DOUBLE PRECISION, INTENT(IN), DIMENSION(m) :: x
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: y
      DOUBLE PRECISION tmp
      DO i = 1,n
        tmp=0.0d0
        DO j = 1,m
          tmp=tmp+a(i,j)*x(j)
        ENDDO
        y(i)=tmp
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE maxv1(a,x,n,m)
c      mat vec mult  x=A*x (sol. stored in input vec)
      IMPLICIT NONE
      INTEGER n,m,i,j
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n,m) :: a
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(m) :: x
      DOUBLE PRECISION,  DIMENSION(m) :: tmp_vec 
      DOUBLE PRECISION tmp
      call vmove(x,tmp_vec,m)
      call vzero(x,n)
      DO i = 1,n
        tmp=0.0d0
        DO j = 1,m
          tmp=tmp+a(i,j)*tmp_vec(j)
        ENDDO
        x(i)=tmp
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE maxb(a,b,c,na,nb,nc,n1,n2,n3)
c      mat mat mult  C=A*B 
      IMPLICIT NONE
      INTEGER na,nb,nc,n1,n2,n3,i,j,k
      DOUBLE PRECISION, INTENT(IN), DIMENSION(na,*) :: a
      DOUBLE PRECISION, INTENT(IN), DIMENSION(nb,*) :: b
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nc,*) :: c
      DOUBLE PRECISION  tmp
      DO i=1,n1
        DO j=1,n2
          tmp = 0.0d0
          DO k=1,n3
            tmp = tmp+a(i,k)*b(k,j)
          ENDDO
          c(i,j)=tmp
        ENDDO
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE matxb(a,b,c,na,nb,nc,n1,n2,n3)
c      mat mat mult  C=At*B (A is transposed)
      IMPLICIT NONE
      INTEGER na,nb,nc,n1,n2,n3,i,j,k
      DOUBLE PRECISION, INTENT(IN), DIMENSION(na,*) :: a
      DOUBLE PRECISION, INTENT(IN), DIMENSION(nb,*) :: b
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nc,*) :: c
      DOUBLE PRECISION tmp
      DO i=1,n1
        DO j=1,n2
          tmp = 0.0d0
          DO k=1,n3
            tmp=tmp+a(k,i)*b(k,j)
          ENDDO
          c(i,j)=tmp
        ENDDO
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE vamb(a,b,c,n)
c      vector operation c = a - b
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: a
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: b
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: c
      DO i=1,n
        c(i)=a(i)-b(i)
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE vamb1(a,b,n)
c      vector operation a = a - b (sol. stored into a)
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: a
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: b
      DOUBLE PRECISION tmp
      DO i=1,n
        tmp=a(i)-b(i)
        a(i)=tmp
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE vmult(a,b,n,c)
c      vec mult by scalar c
c      b = a*c
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: a
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: b
      DOUBLE PRECISION c
      DO i=1,n
        b(i)=a(i)*c
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE vmult1(a,n,c)
c      vec mult by scalar c (sol. stored in input vec)
c      a = a*c
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: a
      DOUBLE PRECISION c,tmp
      DO i=1,n
        tmp=a(i)*c
        a(i)=tmp
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE vapb1(a,b,n)
c      vec b added to vec a 
c      a = a + b
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: a
      DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: b
      DOUBLE PRECISION tmp
      DO i=1,n
        tmp=a(i)+b(i)
        a(i)=tmp
      ENDDO
      END SUBROUTINE
c
      SUBROUTINE symeqn(c,nd,n,nc)
c     SUBROUTINE symeqn(c,nd,n,nc,flag)
c      Function is emulated only partialy.
c      C(n,2*n) 
c      in C(1:n,1:n) is stored matrix to be inverted.
c      Inverted matrix is stored in C(1:n,n+1:2*n)
      IMPLICIT NONE
      INTEGER nd,n,nc,t1
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nd,*) :: c

      call inverse0(c(1,1),c(1,n+1),n)
      t1=0
      end subroutine
c
      subroutine matsym(a,nd,n)
c      building symmetric matrix from lower triangular part
      IMPLICIT NONE
      INTEGER nd,n,i,j
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nd,nd) :: a
      DO i=2,n
        DO j=1,i-1
          a(j,i)=a(i,j)
        ENDDO
      ENDDO
      END SUBROUTINE
c
c Linear Algebra START  ---------------------------------------------
      SUBROUTINE inverse0(a,b,n)
      IMPLICIT NONE 
      INTEGER i,j,k,l,m,n,irow
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,n) :: a,b
      DOUBLE PRECISION big, dum 

      DO i = 1, n
        DO j = 1, n
          b(i,j) = 0.0
        END DO
        b(i,i) = 1.0
      END DO

      DO i = 1, n
        ! in case the entry a(i,i) is zero, we need to find a good pivot;
        ! this pivot
        ! is chosen as the largest value on the column i from a(j,i) with
        ! j = 1,n
        big = a(i,i)
        DO j = i,n
          IF (a(j,i).gt.big) THEN
            big = a(j,i)
            irow = j
          END IF
        END DO
        ! interchange lines i with irow for both a() and b() matrices
        IF (big.gt.a(i,i)) THEN
          DO k = 1, n
            ! matrix a()
            dum = a(i,k)
            a(i,k) = a(irow,k)
            a(irow,k) = dum
            ! matrix b()
            dum = b(i,k)
            b(i,k) = b(irow,k)
            b(irow,k) = dum
          END DO
        END IF
        ! divide all entries in line i from a(i,j) by the value a(i,i);
        ! same operation for the identity matrix
        dum = a(i,i)
        DO j = 1,n
          a(i,j) = a(i,j)/dum
          b(i,j) = b(i,j)/dum
        END DO
        ! make zero all entries in the column a(j,i); same operation for
        ! indent()
        DO j = i+1, n
          dum = a(j,i)
          DO k = 1, n
            a(j,k) = a(j,k) - dum*a(i,k)
            b(j,k) = b(j,k) - dum*b(i,k)
          END DO
        END DO
      END DO
      ! substract appropiate multiple of row j from row j-1
      DO i = 1, n-1
        DO j = i+1, n
          dum = a(i,j)
          DO l = 1, n
            a(i,l) = a(i,l)-dum*a(j,l)
            b(i,l) = b(i,l)-dum*b(j,l)
          END DO
        END DO
      END DO
      END SUBROUTINE

      SUBROUTINE inverse(a,c,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !===========================================================
      IMPLICIT NONE 
      INTEGER i,j,k,n
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,n) :: a,c
      DOUBLE PRECISION,  DIMENSION(n,n) :: L,U
      DOUBLE PRECISION,  DIMENSION(n) :: b,d,x
      DOUBLE PRECISION coeff
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
c      L=0.0
c      U=0.0
c      b=0.0
      DO i=1,n
        b(i)=0.d0
        DO j=1,n
          L(i,j)=0.d0
          U(i,j)=0.d0
        ENDDO
      ENDDO

      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do
      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do
      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      END SUBROUTINE inverse
c
      SUBROUTINE maat(a,c,nc,n,con)
      INTEGER nc,n,i,j
      DOUBLE PRECISION  con
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n) :: a
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nc,*) :: c
      DO i=1,n
        DO j=1,i
          c(i,j)=a(i)*a(j)*con
        ENDDO
      ENDDO
      END SUBROUTINE
#endif

c ANSYS subroutine START  | | |
c                         V V V



*deck,userpl           parallel         user                          dc
      subroutine userpl (elem,intpt,mat,ncomp,kfirst,kfsteq,e,nu,dens,
     x prop,d,ktform,timval,timinc,tem,dtem,toffst,flu,dflu,epel,eppl,
     x statev,usvr,epeq,plwork,sigepl,sigrat,depeq,dt)
c    *********************************************************************
c     Tento program byl vytvoren upravou zdrojoveho kodu souboru USERPL.F,
c     ktery je distribuovan spolecne s programem Ansys 8.0.
c
c         *** ansys(r) copyright(c) 2000
c         *** ansys, inc.
c
c  Vstupni promenne:
c     promenna (typ,def,in/out)    popis
c
c     elem     (int,sc,in)         - cislo prvku
c     intpt    (int,sc,in)         - cislo integracniho bodu daneho prvku
c     mat      (int,sc,in)         - referencni cislo materialu
c     ncomp    (int,sc,in)         - pocet slozek tenzoru napeti/deformace
c                                    4 - x,y,z,xy
c                                    6 - x,y,z,xy,yz,xz
c     kfirst   (int,sc,in)         - 1 pri prvnim pouziti kodu, 0 jindy
c                                    (lze vyuzit pro inicializaci stavovych
c                                     promennych)
c     kfsteq   (int,sc,in)         - 1 pro prvni iteraci daneho kroku NR
c                                    metody, 0 jindy
c
c     e        (dp,sc,in)          - modul pruznosti v tahu                           	1
c     nu       (dp,sc,in)          - poissonovo cislo                                 	1
c     dens     (dp,sc,in)          - hustota materialu                                	0    
c     prop                         - pole obsahujici linearni vlastnosti mat. 		0
c              (dp,ar(9),in)         (ex,ey,ez, gxy,gyz,gxz, nuxy,nuyz,nuxz)			0
c              (dp,ar(1),in)         kdyz ncomp=1 (ex)					0
c     d        (dp,ar(ncomp,ncomp),in)- materialova matice tuhosti				0 sestavi se uvnitr
c     ktform   (int,sc,in)         - klic indikujici vyzadani elastoplasticke		0
c                                    matice
c                                    (=1, vytvorit  .ne.1, nevytvaret)
c
c     timval   (dp,sc,in)          - aktualni hodnota casu					0
c     timinc   (dp,sc,in)          - casovy prirustek v danem kroku NR met.			0
c
c     tem      (dp,sc,in)          - teplota na konci kroku					0
c     dtem     (dp,sc,in)          - prirustek teploty behem kroku				0
c     toffst   (dp,sc,in)          - odchylka teploty od absolutni nuly			0
c     flu      (dp,sc,in)          - prutok na konci kroku					0
c     dflu     (dp,sc,in)          - prirustek prutoku behem kroku				0
c	
c     epel     (dp,ar(ncomp),inout)- modifikovana celkova deformace (testovaci) 		1
c                                    epel = eptot - eppl - eptherm - ...
c     eppl     (dp,ar(ncomp),inout)- plasticka deformace v predchozim kroku		        1
c
c     h   (dp,ar(ncomp,6),inout)- stavove promenne z predchoziho kroku			1
c     usvr     (dp,ar(nuval,nintp),inout)- pridavne stavove promenne z predchozi            ?
c                                    iterace v rovnicich rovnovahy (pro pouziti
c                                    je nutno nadefinovat prikazem nsvr)
c
c     epeq     (dp,sc,inout)       - ekvivalentni plasticka deformace z predch.             1
c                                    kroku
c     plwork   (dp,sc,inout)       - akumulovana plasticka prace z predchoziho		1
c                                    kroku
c
c  Vystupni promenne:
c     promenna (typ,def,in/out)    popis
c
c     epel     (dp,ar(ncomp),inout)- elasticka deformace                                    1
c     eppl     (dp,ar(ncomp),inout)- aktualizovany vektor plasticke deformace               1
c
c     statev   (dp,ar(ncomp,6),inout)- aktualizovane stavove promenne                       1
c     usvr     (dp,ar(nuval,nintp),inout)- aktual. pridavne stavove promenne                1
c
c     epeq     (dp,sc,inout)       - aktualizovana ekvivalentni plasticka def.              1
c     plwork   (dp,sc,inout)       - aktualizovana plasticka prace				1
c
c     sigepl   (dp,sc,out)         - napeti v tahovem diagramu odpovidajici epeq        	1
c     sigrat   (dp,sc,out)         - pomer testovaciho napeti ku mezi kluzu			1
c     depeq    (dp,sc,out)         - prirustek ekvivalentni plasticke deformace		1
c                                    (pouzivano pro auto time stepping - casovy
c                                    krok je redukovan, kdyz depeq dosahne .05)
c
c     dt       (dp,ar(ncomp,ncomp),out)- elastoplasticka matice (operator)			1
c
c  Dalsi parametery  (souvisejici s pouzitim uzivatelskych voleb v Ansysu):
c     promenna (typ,def,in/out)    popis
c     numinp   (int)               - pocet polozek v tabulce materialovych dat
c                                    (prikaz tbdat)
c     nuval    (int)               - pocet pridavnych stavovych promennych
c                                    na jeden integracni bod
c     nintp    (int)               - maximalni pocet integracnich bodu daneho
c                                    prvku pro pouziti s danou subrutinou
c                                    (14 je maximum)
c        pozor:  nuval x nintp = nstv(v prikazu nsvr); nesmi presahnout 840!
c
c  Vnitrni promenne:
c     promenna (typ,def,in/out)  popis
c     b        (dp,ar(6,6))    - pomocna matice
c     c        (dp,ar(6,12))   - pomocna matice (take pro inverzi matice)
c     con      (dp,sc)         - pomocna promenna
c     deppl    (dp,ar(6))      - prirustek vektoru plasticke deformace
c     dfds     (dp,ar(6))      - derivace funkce plasticity (smerovy vektor n)
c     dlamb    (dp,sc)         - plasticky multiplikator
c     h        (dp,sc)         - tecny plasticky modul (napeti/plasticka deformace)
c     n2       (int,sc)        - ncomp*ncomp
c     seqtr    (dp,sc)         - eqvivalentni (von mises) testovaci napeti
c     sigtr    (dp,ar(6))      - testovaci napeti
c     sigy     (dp,sc)         - mez kluzu
c     vect     (dp,ar(6))      - pomocny vektor
c    *************************pridane promenne*********************
c     promM  (int)                    - pocet kinematickych casti (backstress)
c     k      (int)                    - pocitadlo
c     i      (int)                    - pocitadlo
c     j      (int)                    - pocitadlo
c     con2   (dp,sc)                  - pomocna promenna
c     g      (dp,sc)                  - modul pruznosti ve smyku
c     ni     (dp,sc)                  - mat. konstanta ni1nek pro ratchetting
c     mrat   (dp,sc)                  - mat. konstanta pro multiaxialni ratchetting
c     pk     (dp,sc)                  - akumulovana pl. def. v k-te iteraci
c     dlam0  (dp,sc)                  - dlambda v (k-1)-ve iteraci
c     dlam1  (dp,sc)                  - dlambda v (k-2)-ve iteraci
c     a      (dp,ar(ncomp,promM))     - matice obsahujici slozky kin. tenzoru
c     an     (dp,ar(6))               - backstress v n+1-tem kroku
c     ano    (dp,ar(6))               - backstress v n-tem kroku
c     sa     (dp,ar(6))               - rozdil s-a v kroku n+1
c     r      (dp,ar(promM))           - mat. konstanty
c     gamai  (dp,ar(promM))           - mat. konstanty
c     ci     (dp,ar(promM))           - mat. konstanty
c     fi     (dp,ar(promM))           - hledana promenna v iteracni procedure
c     hf     (int,ar(6))              - identifikator
c     hi     (dp,ar(6,6))             - matice Hn+1 pro i-ty backstress
c     b      (dp,ar(6,6))             - pomocna matice
c     bi     (dp,ar(6,6))             - pomocna matice
c     c      (dp,ar(6,6))             - pomocna matice
c     eps2   (dp,sc)                  - kriterium pro Aitkenovu formuli
c     conb   (dp,sc)                  - konstanta b
c     Rnek   (dp,sc)                  - konstanta Rnek
c     omega  (dp,sc)                  - konstanta omega
c     p      (dp,sc)                  - konstanta p....akumulovana pl.def.
c     conk   (dp,sc)                  - konstanta k ... odpovida sigY na zacatku
c     
#ifdef inAnsys
#include "impcom.inc"
      external erhandler
#endif

c
c  userd defined parameters
c     --- viz prikaz Ansysu nsvr
c     --- see Ansys command 'nsvr'
      integer nuval,nintp,iott
      parameter (nuval=1,nintp=1)
c     --- number of entities in table of definition of material 
c                                             (command tbdat)
      integer numinp,promM
      parameter (numinp=22,promM=6)
      double precision eps,eps2
      parameter (eps=1.0d-8,eps2=1.0d-4)
c
#ifdef inAnsys
c  exter subroutines and functions
      external vdot
      external tbuser,vmove,vamb,maxv,vmult,vamb1,vapb1,
     x vmult1,vzero,maat,symeqn,matsym,maxb,maxv1,matxb
#endif
      double precision vdot
c
c   integers 
      integer elem,intpt,mat,ncomp,kfirst,kfsteq,ktform,
     x n2,k,i,j,hf(promM)
c
c   double precision
      double precision e,nu,g,dens,prop(*),d(ncomp,ncomp),timval,timinc,
     x tem,dtem,toffst,flu,dflu,epel(ncomp),eppl(ncomp),statev(ncomp,6),
     x epeq,plwork,sigepl,sigrat,depeq,b(ncomp,ncomp),bi(ncomp,ncomp),
     x dt(ncomp,ncomp),proptb(numinp),con,deppl(ncomp),dfds(ncomp),
     x ano(ncomp),seqtr,sigtr(ncomp),sigy,vect(ncomp),pk,a(ncomp,promM),
     x r(promM),gamai(promM),ci(promM),fi(promM),ni,an(ncomp),con2,h,
     x dlam1,c(6,12),dlamb,dlam0,hi(ncomp,ncomp),sa(ncomp),
     x usvr(nuval,nintp),conb1,p,ni0,nas(promM),mrat,
     x nis(promM),ni02,ni2,conb2,conb,Rnek,conk
c
#ifdef inAnsys
c  due to FORTRAN: Conditional tracking for qa activated by stability
#ifdef TRACK_QA
      external TrackPoint
      call TrackPoint('userpl')
#endif
c
      if (nuval*nintp.gt.840)
     X call erhandler('userpl',5010,4,
     x 'Maximum storage allowed by USVR has been exceeded.'
     x ,0.0d0,' ')
#endif
c
c --- inicializace elastoplasticke matice pro pripad elasticity
      n2 = ncomp*ncomp
      call vmove (d(1,1),dt(1,1),n2)
      depeq = 0.0d0
      g = 0.5d0*e / (1 + nu)
c *** pomocne matice
      call vzero (b(1,1),n2)
      do j=1,3
        b(j,j) = 1.0d0
      end do
      do j=4,ncomp
        b(j,j) = 2.0d0
      end do
c *** reading of new variables
      call vzero (an(1),ncomp)
      do i=1,promM
        fi(i) = 1.0d0
        do j=1,ncomp
          a(j,i) = statev(j,i)
          an(j) = an(j) + a(j,i)
        end do
      end do
c *** getting user defined material parameters
      call tbuser (mat,numinp,proptb(1))
      conk = proptb(1)
      ni0 = proptb(2)
      conb1 = proptb(3)
      ni = proptb(4)
      do i=1,promM
        ci(i) = proptb(2*i+3)
        gamai(i) = proptb(2*i+4)
        r(i) = ci(i)/gamai(i)
      end do
      mrat = proptb(17)
      ni02 = proptb(18)
      conb2 = proptb(19)
      ni2 = proptb(20)
      conb = proptb(21)
      Rnek = proptb(22)
      sigy = conk + Rnek*(1.0d0-exp(-conb*epeq))
c *** another control of correctness of the calculation of constant  
c choice added
      if (sigy.le.0.0d0) go to 999
      if (r(6).le.0.0d0) go to 999
      iott = 6
c
c *** testing stress calculation
c
      call vmove (an(1),ano(1),ncomp)
      call maxv (d(1,1),epel(1),sigtr(1),ncomp,ncomp)
      dlam0 = 0.0d0
      dlam1 = 0.0d0
      p = epeq
      con = (sigtr(1)+sigtr(2)+sigtr(3))/3.0d0
      sigtr(1) = sigtr(1) - con
      sigtr(2) = sigtr(2) - con
      sigtr(3) = sigtr(3) - con
      call vamb (sigtr(1),ano(1),vect(1),ncomp)
c *** von Mises stress calculation
      seqtr = vect(1)**2 + vect(2)**2 + vect(3)**2 + 2.0d0*vect(4)**2
      if (ncomp.eq.6) seqtr = seqtr + 2.0d0*(vect(5)**2 + vect(6)**2)
      seqtr = 1.5d0*seqtr
c *** ratio to yield stress
      sigrat = sqrt(seqtr/(sigy*sigy))
c
c *** plastic state checking
c
      con = sigrat - 1.0d0
      if (con.lt.eps) then
        go to 999
      endif
c *** loop for succesive substitution
      do k=1,500
        call vzero (an(1),ncomp)
        do j=1,ncomp
          do i=1,promM
            an(j) = an(j) + fi(i)*a(j,i)
          end do
        end do
        do j=1,ncomp
          sa(j) = sigtr(j) - an(j)
        end do        
        con2 = 0.0d0
        do i=1,promM
          con2 = con2 + fi(i) * ci(i)
        end do
        con2 = con2 + 3.0d0 * g
        con = sa(1)**2+sa(2)**2+sa(3)**2+2.0d0*sa(4)**2
        if (ncomp.eq.6) con = con + 2.0d0*(sa(5)**2 + sa(6)**2)
        con = sqrt(1.5d0*con)
        dlamb = (con - sigy)/con2
c *** Aitkenov's formula
        if (MOD(k,3).eq.0) then
          con=dlamb-(dlamb-dlam0)*(dlamb-dlam0)/(dlamb-2.0d0*dlam0+
     x          dlam1)
          if (con.gt.0.0d0) dlamb=con
        endif
        con = sigy/(sigy + dlamb*con2)
        call vmult1 (sa(1),ncomp,con)
c *** calculation of increment of plastic deformation
        con = 1.5d0*dlamb/sigy
        call vmult (sa(1),deppl(1),ncomp,con)
        call maxv1 (b(1,1),deppl(1),ncomp,ncomp)
        con2 = abs(dlamb-dlam0)/dlamb
        p = epeq + dlamb
        con = sqrt(1.5d0)/sigy
        call vmult (sa(1),dfds(1),ncomp,con)
        sigy = conk + Rnek*(1.0d0-exp(-conb*p))
        if (con2.lt.eps2) EXIT
        dlam1 = ni + (ni0 - ni)*exp(-conb1*p) +
     x            ni2 +(ni02 - ni2)*exp(-conb2*p)
        do i=1,promM
          do j=1,3
            vect(j) = a(j,i)+2.0d0*ci(i)*deppl(j)/3.0d0
          end do
          do j=4,ncomp
            vect(j) = a(j,i)+ci(i)*deppl(j)/3.0d0
          end do
          con = vect(1)**2+vect(2)**2+vect(3)**2+2.0d0*vect(4)**2
          if (ncomp.eq.6) con=con+2.0d0*(vect(5)**2 + vect(6)**2)
          con = sqrt(1.5d0*con)
          if (con.ne.0.0d0) then
            call vmult1 (vect(1),ncomp,1.0d0/con)
            nas(i) = 0.0d0
            do j=1,ncomp
              nas(i) = nas(i)+sqrt(1.5d0)*dfds(j)*vect(j)
            end do 
            nas(i) = abs(nas(i))
          else
            nas(i) = 1.0d0
          endif
          nas(i)=nas(i)**mrat
          nis(i) = dlam1*nas(i)
          con2 = r(i)/(r(i)+nis(i)*ci(i)*dlamb)
          do j=1,3
            vect(j) = con2*(a(j,i)+(2.0d0*ci(i)*deppl(j)/3.0d0))
          end do
          do j=4,ncomp
            vect(j) = con2*(a(j,i)+(ci(i)*deppl(j)/3.0d0))
          end do
          con = vect(1)**2+vect(2)**2+vect(3)**2+2.0d0*vect(4)**2
          if (ncomp.eq.6) con=con+2.0d0*(vect(5)**2 + vect(6)**2)
          con = sqrt(1.5d0*con)
          if (con.gt.r(i)) then
            fi(i) = r(i)*con2/con
            hf(i) = 1
          else
            fi(i) = con2
            hf(i) = 0
          endif
        end do
        dlam1 = dlam0
        dlam0 = dlamb
      end do
      call vzero (an(1),ncomp)
c *** state variables storing
      do i=1,promM
        do j=1,3
          con2 = fi(i)*(a(j,i)+(2.0d0*ci(i)*deppl(j)/3.0d0))
          a(j,i) = con2
          statev(j,i) = con2
          an(j) = an(j) + con2
        end do
        do j=4,ncomp
          con2 = fi(i)*(a(j,i)+(ci(i)*deppl(j)/3.0d0))
          a(j,i) = con2
          statev(j,i) = con2
          an(j) = an(j) + con2
        end do
      end do
c *** update of deformations and other variables
c
      call vapb1 (eppl(1),deppl(1),ncomp)
      call vamb1 (epel(1),deppl(1),ncomp)
      depeq = dlamb
      epeq = epeq + depeq
      call vamb (an(1),ano(1),sigtr(1),ncomp)
      call vmult1 (sigtr(1),ncomp,0.5d0)
      call maxv (d(1,1),epel(1),vect(1),ncomp,ncomp)
      call vamb (vect(1),sigtr(1),ano(1),ncomp)
      plwork = plwork + vdot(ano(1),deppl(1),ncomp)
      sigepl = (vect(1)-vect(2))**2 + (vect(2)-vect(3))**2 +
     x           (vect(3)-vect(1))**2 + 6.0d0*vect(4)**2
      if (ncomp.eq.6) sigepl = sigepl +
     x      6.0d0*(vect(5)**2 + vect(6)**2)
      sigepl = sqrt (0.5d0*sigepl)
c *** elastoplastic matrix calculation
 500  if (ktform.eq.1) then
c *** consisten matrix calculation
        if (dlamb.gt.1.0d-8) then
          call vzero (hi(1,1),n2)
          call vzero (bi(1,1),n2)
          do j=1,3
            bi(j,j) = 1.0d0
          end do
          do j=4,ncomp
            bi(j,j) = 0.5d0
          end do
          do j=4,ncomp
            b(j,j) = 1.0d0
          end do
          do i=1,promM
            con = sqrt(1.5d0)/r(i)
            con2 = 2.0d0*fi(i)*ci(i)/3.0d0
            do j=1,ncomp
              ano(j) = con*a(j,i)
            end do
            do j=1,ncomp
              do k=1,ncomp
                hi(j,k)=hi(j,k)+con2*(bi(j,k)-nis(i)*ano(j)*dfds(k)
     x            -hf(i)*(ano(j)*ano(k)-nis(i)*ano(j)*dfds(k)))
              end do
            end do                 
          end do
          do j=4,ncomp
            bi(j,j) = 2.0d0
          end do
          con = 2.0d0*sigy/(3.0d0*depeq)
          con2 = (2.0d0/3.0d0)*conb*(Rnek-sigy+conk)
          do j=1,ncomp
            do k=1,ncomp
              c(j,k) = hi(j,k) + (con2-con)*dfds(j)*dfds(k)
            end do
          end do
          call matxb (bi(1,1),c(1,1),hi(1,1),ncomp,6,ncomp,
     x        ncomp,ncomp,ncomp)
          do j=1,ncomp
            do k=1,ncomp
              c(j,k) = hi(j,k) + (2.0d0*g+con)*b(j,k)
            end do
          end do
          call symeqn (c(1,1),6,ncomp,-ncomp)
          call vzero (hi(1,1),n2)
          hi(1,1) = 1.0d0
          hi(2,2) = 1.0d0
          hi(3,3) = 1.0d0
          hi(4,4) = 1.5d0
          if (ncomp.eq.6) then
            hi(5,5) = 1.5d0
            hi(6,6) = 1.5d0
          endif
          hi(2,1) = -0.5d0
          hi(3,2) = -0.5d0
          hi(3,1) = -0.5d0
          con2 = 2.0d0/3.0d0
          call vmult1 (hi(1,1),n2,con2)
          call matsym (hi(1,1),ncomp,ncomp)
          call maxb (c(1,ncomp+1),hi(1,1),b(1,1),6,ncomp,ncomp,
     x              ncomp,ncomp,ncomp)
          call maxb (d(1,1),b(1,1),hi(1,1),ncomp,ncomp,ncomp,
     x         ncomp,ncomp,ncomp)
          con = 2.0d0*g
          do j=1,ncomp
            do k=1,ncomp
              dt(j,k) = d(j,k) - con*hi(j,k)
            end do
          end do
          call matsym (dt(1,1),ncomp,ncomp)
        else
c *** inconsistent modul calculation
          h = 0.0d0
          do i=1,promM
            con = sqrt(1.5d0)/r(i)
            con2 = 0.0d0
            do j=1,ncomp
              ano(j) = con*a(j,i)
              con2 = con2+ano(j)*dfds(j)
            end do
            do j=1,ncomp
              vect(j)=gamai(i)*a(j,i)*(nis(i)+hf(i)*(con2-nis(i)+
     x                abs(con2-nis(i)))/2.0d0)
            end do
            con2 = sqrt(1.5d0)
            con = 0.0d0
            do j=1,ncomp
              con = con + con2*dfds(j)*b(j,j)*vect(j)
            end do
            h = h - con + ci(i)
          end do
          h = h + conb*(Rnek-sigy+conk)
          con = 6.0d0*(g**2)/(3.0d0*g+h)
          call maat (dfds(1),dt(1,1),ncomp,ncomp, -con)
          call matsym (dt(1,1),ncomp,ncomp)
        endif
      endif
c
 999  return
      end
