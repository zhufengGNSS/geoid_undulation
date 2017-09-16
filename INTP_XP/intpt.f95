
    module intptFunc


      contains

      subroutine INTPT( VecFLAT, VecFLON , VecU, len1)
!!!CCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+C-CCCCCCC+
!!!
!!! GRID FILE STRUCTURE:
!!!
!!!  This FORTRAN program will interpolate values from a WORLD-WIDE grid file.
!!!  The grid file must have a header as record 1 with the following format:
!!!
!!!       south    north    west    east    spacing N-S   spacing W-E
!!!
!!!  The grid data are in subsequent records arranged from North to South, West
!!!  to East, i.e. record 2 would be the northern most latitude band of data.
!!!
!!!-------------------------------------------------------------------------------
!!!
!!!  Several interpolation schemes are possible using this program.  By setting
!!!  the parameter IWINDO to a particular value, one may use Bi-linear interpo-
!!!  lation or a spline interpolation with a window of size of IWINDO  x IWINDO.
!!!  Check the subroutine "INTERP" for more information.
!!!
!!!  This program is set up to extend the longitude values by a minimum of 2 posts
!!!  on the positive side of the east longitude limit and a minimum of 2 posts to
!!!  the negative side of the west longitude.
!!!
!!!-------------------------------------------------------------------------------
!!!
!!!  Change the following variables in the parameter statements below to match
!!!  the grid file information and desired interpolation method:
!!!
!!!                IWINDO = Size of window of interpolation and method.
!!!                             IWINDO=0 gives bilinear interpolation.
!!!                             Check subroutine "INTERP" below for more
!!!                             information.
!!!
!!!                ****NOTE: If IWINDO=0 be sure to set NBDR=4 *****
!!!                             In general, NBDR=2*IWINDO
!!!
!!!
!!!                NLAT = Number of values in latitude.
!!!                NLON = Number of values in longitude. (Leave the NBDR
!!!                            parameter in the NLON parameter declaration.
!!!                            This controls the extension of the grid in
!!!                            longitude.)
!!!                DLAT = Spacing in latitude.
!!!                DLON = Spacing in longitude.
!!!
!!!-------------------------------------------------------------------------------
!!!
!!! FILE INFORMATION:
!!!
!!!  INPUT DATA GRID: UNIT = 1
!!!
!!!    For structure see above.
!!!
!!!
!!!  INPUT DATA POINTS:  UNIT = 11
!!!
!!!     The input file of points to be interpolated is arranged as:
!!!
!!!                       latitude   longitude
!!!                         (decimal degrees)
!!!
!!!
!!! OUPUT OF INTERPOLATED POINTS:  UNIT = 20
!!!
!!!    The output file of interpolated points is arranged as:
!!!
!!!                       latitude   longitude    value
!!!                        (decimal degrees)
!!!
!!!-------------------------------------------------------------------------------
!!!
!!!  This program was provided by Professor Richard H. Rapp of The Ohio State
!!!    University, December, 1996.  It was put into its present form by the National
!!!    Imagery and Mapping Agency (NIMA), December, 1996.
!!!
!!!CCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+C-CCCCCCC+

      use ISO_C_BINDING

      implicit real*8(a-h,o-z)
      real*8 south,north,west,east,dphi,dlam
      integer(c_int) :: len1
      real(c_double),  intent(in) :: VecFLAT(len1), VecFLON(len1)
      real(c_double),  intent(out) :: VecU(len1)


      PARAMETER(IWINDO=4)
      PARAMETER(NBDR=2*IWINDO)
      PARAMETER(NLAT=721,NLON=(1441+NBDR),DLAT=0.25D0,DLON=0.25D0)
      DIMENSION A1(NLON)
      REAL*4 H (nlat, nlon)
      CHARACTER NAM*41
      REAL*8 XMIN,XMAX,RDX
      DATA XMIN, XMAX, RDX, NAM/-50.,50., 2.,&
           'GEOID HTS'/

!!!      CALL HISTO(6,1,NAM,XMIN,XMAX,RDX,1,0.D0,1)

!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCC                                                           CCCCCCCCCCC
!!!CCCCCCCCC                 INPUT/OUPUT FILE NAMES                    CCCCCCCCCCC
!!!CCCCCCCCC            MAKE CHANGES TO FILE NAMES BELOW               CCCCCCCCCCC
!!!CCCCCCCCC                                                           CCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!!! data grid file
      open(1,file='WW15MGH.GRD',status='old')

!!! input file of points to be interpolated
      open(11,file='INPUT.DAT',status='old')

!!! output file of interpolated values
      open(20,file='OUTINTPT.DAT')


!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!!! read header of input grid file.
      read(1,*)south,north,west,east,dphi,dlam
!!!     write(*,*)'input file header:'
!!!     write(*,998)south,north,west,east,dphi,dlam,
!!!    &  int((north-south)/dphi+1),int((east-west)/dlam+1)
!!!     write(*,*)'Southwest point in grid storage array:'
!!!     write(*,998)south,west-nbdr/2*dlam
998   format(6f8.2,2i6)


!!! read input grid file and store in array H
      DO I=1,NLAT
         READ(1,*) (A1(K),K=1,nlon-nbdr)
         DO J=1,nlon-nbdr
            H(NLAT-I+1,J+nbdr/2)=A1(J)
         ENDDO
         DO K=1,nbdr/2
            H(NLAT-I+1,K)=A1(nlon-nbdr-nbdr/2+K)
3           H(NLAT-I+1,k+nlon-nbdr/2)=A1(K)
         ENDDO
      ENDDO

!!! read in the points to be interpolated
 onePoint :   do Nindex = 1, size(VecFLAT, 1)
        FLAT = VecFLAT(Nindex)
        FLON = VecFLON(Nindex)
!!! 80    READ(11,*,END=100) FLAT,FLON

!!! call interpolation subroutine
      CALL INTERP(IWINDO,12.D0,H,south,west-nbdr/2*dlam,DLAT,DLON,&
                  NLAT,NLON,NLAT,NLON,FLAT,FLON,UN)

!!!      CALL HISTO(6,2,NAM,XMIN,XMAX,RDX,1,UN,1)

!!! output the interpolated points
    !!!    WRITE(20,997)FLAT,FLON,UN
    VecU(Nindex) = UN
997   format(2f14.7,f12.3)

      !!! GOTO 80
    end do onePoint

100   CONTINUE

!!!      CALL HISTO(6,3,NAM,XMIN,XMAX,RDX,1,0.D0,1)

!!!  end of the main program
      return
      END subroutine



      SUBROUTINE INTERP(IWINDO,DMIN,H,PHIS,DLAW,DDFI,DDLA,NPHI,NDLA,&
                        IPDIM,ILDIM,PHI,DLA,VALINT)
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!                                                                      C
!!!     SUBROUTINE FOR INTERPOLATION OF VALUES FROM A STANDARD DTM-GRID  C
!!!     TO INDIVIDUAL STATION LOCATIONS.                                 C
!!!                                                                      C
!!!                                                                      C
!!!     INPUT PARAMETERS...                                              C
!!!     ===================                                              C
!!!     IWINDO...    A SPLINE WINDOW OF SIZE 'IWINDO' X 'IWINDO' WILL BE C
!!!                  USED AROUND EACH STATION. IF 'IWINDO' IS 0 OR 1,    C
!!!                  BILINEAR INTERPOLATION WILL BE USED.                C
!!!     DMIN...      MINIMUM ACCEPTABLE DISTANCE FROM THE GRID EDGE IN   C
!!!                  KM (USEFUL FOR FFT GRIDS).                          C
!!!     H...         2D DATA ARRAY (ELEMENT (1,1) IN SW CORNER).         C
!!!     PHIS,DLAW... LATITUDE AND LONGITUDE OF SW GRID POINT.            C
!!!     DDFI,DDLA... GRID SPACING IN LATITUDE AND LONGITUDE DIRECTION.   C
!!!     NPHI,NDLA... NUMBER OF GRID POINTS IN LATITUDE AND LONGITUDE     C
!!!                  DIRECTION.                                          C
!!!     IPDIM,ILDIM..DIMENSIONS OF 2D DATA ARRAY 'H' AS DECLARED IN THE  C
!!!                  CALLING PROGRAM.                                    C
!!!     PHI,DLA...   LATITUDE AND LONGITUDE OF INTERPOLATION POINT.      C
!!!                                                                      C
!!!                                                                      C
!!!     OUTPUT PARAMETERS...                                             C
!!!     ====================                                             C
!!!     VALINT...    INTERPOLATED VALUE.                                 C
!!!                                                                      C
!!!                                                                      C
!!!     EXECUTION TIME ON CDC 990 IS...                                  C
!!!     ===============================                                  C
!!!     +------------------+-------------------+-------------------+     C
!!!     I  INTERPOLATION   I  OPT=LOW          I  OPT=HIGH         I     C
!!!     I------------------I-------------------I-------------------I     C
!!!     I  BILINEAR        I  1.44 MSEC/STAT.  I  1.44 MSEC/STAT.  I     C
!!!     I  3 X 3 SPLINE    I  1.53 MSEC/STAT.  I  1.51 MSEC/STAT.  I     C
!!!     I  5 X 5 SPLINE    I  1.70 MSEC/STAT.  I  1.67 MSEC/STAT.  I     C
!!!     I  7 X 7 SPLINE    I  2.02 MSEC/STAT.  I  1.74 MSEC/STAT.  I     C
!!!     I  9 X 9 SPLINE    I  2.31 MSEC/STAT.  I  2.00 MSEC/STAT.  I     C
!!!     +------------------+-------------------+-------------------+     C
!!!                                                                      C
!!!                                                                      C
!!!     PROGRAM CREATION BY...   H. DENKER          MAY 30, 1987         C
!!!                              H. DENKER          MARCH 13, 1989       C
!!!                                                                      C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER(IPA1=20)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LODD
      REAL*4 H(IPDIM,ILDIM)
      DIMENSION A(IPA1),R(IPA1),Q(IPA1),HC(IPA1)
      IDIM1=IPA1
      RHO=57.29577951D0
      REARTH=6371000.D0
      IF(IWINDO.LT.2) IWINDO=2
      IF(IWINDO.GT.IDIM1) IWINDO=IDIM1
      ILIM=DMIN*1000.*RHO/(REARTH*DDFI)
      JLIM=DMIN*1000.*RHO/(REARTH*DDLA*COS((PHIS+DDFI*NPHI/2.)/RHO))
      LODD=(IWINDO/2)*2.NE.IWINDO
      RI=(PHI-PHIS)/DDFI
      RJ=(DLA-DLAW)/DDLA
      IF(LODD) THEN
        I0=RI-0.5
        J0=RJ-0.5
      ELSE
        I0=RI
        J0=RJ
      ENDIF
      I0=I0-IWINDO/2+1
      J0=J0-IWINDO/2+1
      II=I0+IWINDO-1
      JJ=J0+IWINDO-1
      IF(I0.LT.0 .OR. II.GE.NPHI .OR. J0.LT.0 .OR. JJ.GE.NDLA) THEN
        WRITE(6,7008) PHI,DLA
        VALINT=999999.
        RETURN
      ELSEIF(I0.LT.ILIM .OR. II.GT.NPHI-ILIM .OR. J0.LT.JLIM .OR.&
        JJ.GT.NDLA-JLIM) THEN
!!!        IF(NPOINT.LE.ILIST) WRITE(6,7009) PHI,DLA
        VALINT=999999.
        RETURN
      ENDIF
7008  FORMAT(' ',2F10.6,' STATION TOO NEAR GRID BOUNDARY  - NO INT.'&
      ,' POSSIBLE|')
7009  FORMAT(' ',2F10.6,' STATION OUTSIDE ACCEPTABLE AREA - NO INT.'&
      ,' PERFORMED|')
      IF(IWINDO.GT.2) THEN
        DO 110 I=1,IWINDO
          DO 111 J=1,IWINDO
            A(J)=H(I0+I,J0+J)
111       CONTINUE
          CALL INITSP(A,IWINDO,R,Q)
          HC(I)=SPLINE(RJ-J0+1.,A,IWINDO,R)
110     CONTINUE
        CALL INITSP(HC,IWINDO,R,Q)
        VALINT=SPLINE(RI-I0+1.,HC,IWINDO,R)
      ELSE
        VALINT=BILIN(RI+1.,RJ+1.,H,NPHI,NDLA,IPDIM,ILDIM)
      ENDIF
      RETURN
      END subroutine



      FUNCTION BILIN(RI,RJ,A,IMAX,JMAX,IADIM,JADIM)
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!                                                                      C
!!!                           B I L I N                                  C
!!!                                                                      C
!!!  INTERPOLATES VALUES IN AN ARRAY A USING BILINEAR                    C
!!!  (PARABOLIC HYPERBOLOID) INTERPOLATION.                              C
!!!                                                                      C
!!!----------------------------------------------------------------------C
!!!                                                                      C
!!!  PARAMETERS:                                                         C
!!!                                                                      C
!!!  BILIN...       INTERPOLATED VALUE                                   C
!!!                                                                      C
!!!  RI, RJ...      INTERPOLATION ARGUMENT, (1,1) IN LOWER LEFT CORNER,  C
!!!                 (IMAX, JMAX) IN UPPER RIGHT.                         C
!!!                                                                      C
!!!  A...           INTEGER*2 ARRAY WITH ARGUMENTS                       C
!!!                                                                      C
!!!  IMAX, JMAX...  NUMBER OF POINTS IN GRID                             C
!!!                                                                      C
!!!  IADIM, JADIM...DECLARED DIMENSIONS OF 'A'                           C
!!!                                                                      C
!!!  OUTSIDE AREA COVERED BY 'A' THE FUNCTION RETURNS THE VALUE OF       C
!!!  THE NEAREST BOUNDARY POINT.                                         C
!!!                                                                      C
!!!----------------------------------------------------------------------C
!!!                                                                      C
!!!  PROGRAMMER:                                                         C
!!!  RENE FORSBERG, JULY 1983                                            C
!!!                                                                      C
!!!  MODIFICATIONS BY:                                                   C
!!!  HEINER DENKER, 07/01/1987                                           C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION A(IADIM, JADIM)
      REAL*4 A
      IN = IFRAC(RI)
      IE = IFRAC(RJ)
      RN = RI - IN
      RE = RJ - IE
      IF (IN.LT.1) THEN
        IN = 1
        RN = 0.0
      ELSEIF (IN.GE.IMAX) THEN
        IN = IMAX-1
        RN = 1.0
      ENDIF
      IF (IE.LT.1) THEN
        IE = 1
        RE = 0.0
      ELSEIF (IE.GE.JMAX) THEN
        IE = JMAX-1
        RE = 1.0
      ENDIF
      RNM1=1.-RN
      REM1=1.-RE
      BILIN = RNM1*REM1*A(IN,IE) +&
      RN*REM1*A(IN+1,IE) + RNM1*RE*A(IN,IE+1) +&
      RN*RE*A(IN+1,IE+1)
      RETURN
      END FUNCTION


      SUBROUTINE INITSP(Y, N, R, Q)
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!                                                                      C
!!!                      I N I T S P                                     C
!!!                                                                      C
!!!  INITIALIZATION PROCEDURE FOR FAST 1-DIMENSIONAL EQUIDISTANT         C
!!!  SPLINE INTERPOLATION, WITH FREE BOUNDARY END CONDITIONS             C
!!!  REFERENCE: JOSEF STOER: EINFUHRUNG IN DIE NUMERISCHE MATHEMATIK     C
!!!  I, SPRINGER 1972, PAGE 82 AND 86.                                   C
!!!                                                                      C
!!!----------------------------------------------------------------------C
!!!                                                                      C
!!!  PARAMETERS (REAL):                                                  C
!!!                                                                      C
!!!  Y...   GIVEN VALUES, Y(1), ..., Y(N)                                C
!!!                                                                      C
!!!  R...   SPLINE MOMENTS (1 ... N), TO BE USED BY FUNCTION 'SPLINE'    C
!!!                                                                      C
!!!  Q...   WORK-ARRAY, DECLARED AT LEAST 1:N                            C
!!!                                                                      C
!!!----------------------------------------------------------------------C
!!!                                                                      C
!!!  RENE FORSBERG, JULY 1983                                            C
!!!                                                                      C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N), R(N), Q(N)
      Q(1) = 0.0
      R(1) = 0.0
      DO 11 K = 2, N-1
        P = Q(K-1)/2+2
        Q(K) = -0.5/P
        R(K) = (3*(Y(K+1)-2*Y(K)+Y(K-1)) - R(K-1)/2)/P
   11 CONTINUE
      R(N) = 0.0
      DO 12 K = N-1, 2, -1
        R(K) = Q(K)*R(K+1)+R(K)
   12 CONTINUE
      RETURN
      END subroutine




      FUNCTION SPLINE(X, Y, N, R)
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!                                                                      C
!!!                          S P L I N E                                 C
!!!                                                                      C
!!!  FAST ONE-DIMENSIONAL EQUIDISTANT SPLINE INTERPOLATION FUNCTION.     C
!!!  REFERENCE: JOSEF STOER: EINFUHRUNG IN DIE NUMERISCHE MATHEMATIK     C
!!!  I, SPRINGER 1972, PAGE 81.                                          C
!!!                                                                      C
!!!----------------------------------------------------------------------C
!!!                                                                      C
!!!  PARAMETERS:                                                         C
!!!                                                                      C
!!!  X...  INTERPOLATION ARGUMENT (REAL), X = 1 FIRST DATA-POINT,        C
!!!        X = N LAST DATA-POINT. OUTSIDE THE RANGE LINEAR EXTRA-        C
!!!        POLATION IS USED.                                             C
!!!                                                                      C
!!!  Y...  REAL*8 ARRAY, 1 .. N : DATA VALUES                            C
!!!                                                                      C
!!!  R...  DO: SPLINE MOMENTS CALCULATED BY SUBROUTINE 'INITSP'          C
!!!                                                                      C
!!!----------------------------------------------------------------------C
!!!                                                                      C
!!!  PROGRAMMER:                                                         C
!!!  RENE FORSBERG, JUNE 1983                                            C
!!!                                                                      C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H, P-Z)
      DIMENSION Y(N), R(N)
      IF (X.LT.1) THEN
        SPLINE = Y(1) + (X-1)*(Y(2)-Y(1)-R(2)/6)
      ELSEIF (X.GT.N) THEN
        SPLINE = Y(N) + (X-N)*(Y(N)-Y(N-1)+R(N-1)/6)
      ELSE
        J = IFRAC(X)
        XX = X - J
        SPLINE = Y(J) +&
                 XX * ((Y(J+1)-Y(J)-R(J)/3-R(J+1)/6) +&
                 XX * (R(J)/2 +&
                 XX * (R(J+1)-R(J))/6))
      ENDIF
      RETURN
      END FUNCTION


      FUNCTION IFRAC(R)
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!                                                                      C
!!!                    F U N C T I O N   I F R A C                       C
!!!                  ===============================                     C
!!!                                                                      C
!!!  SUBROUTINE GIVING TRUE INTEGER PART OF A REAL E.G.                  C
!!!                                                                      C
!!!    FOR   1. = R < 2.   IS    IFRAC = 1                               C
!!!    FOR   0. = R < 1.   IS    IFRAC = 0                               C
!!!    FOR  -1. = R < 0.   IS    IFRAC =-1                               C
!!!    FOR  -2. = R <-1.   IS    IFRAC =-2                               C
!!!                                                                      C
!!!  RF, JUNE 1983                                                       C
!!!  HD, JANUARY 1987                                                    C
!!!                                                                      C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8(A-H,O-Z)
      IFRAC=R
      IF (R.GE.0) RETURN
      IF (R.EQ.IFRAC) RETURN
      IFRAC = IFRAC - 1
      RETURN
      END FUNCTION


      SUBROUTINE HISTO(IUN6,IENT,NAM,XMIN,XMAX,RDX,J,X,NEWP)
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!                                                                      C
!!!     SUBROUTINE HISTO, VERSION 830909 FTN5.                           C
!!!                                                                      C
!!!     THE SUBROUTINE HISTO COMPUTES UP TO 20 HISTOGRAMS FOR DIFFERENT  C
!!!     PARAMETERS IN PARALLEL EXECUTION AND PRINTS THE HISTOGRAMS       C
!!!     ON LINE PRINTER.                                                 C
!!!                                                                      C
!!!     INPUT PARAMETER DESCRIPTION...                                   C
!!!     ==============================                                   C
!!!                                                                      C
!!!     IUN6...      UNIT NUMBER OF LINE PRINTER OUTPUT. THIS UNIT HAS   C
!!!                  TO BE DECLARED AS OUTPUT IN THE PROGRAM CARD.       C
!!!                  (EXAMPLE...IUN6=6, DECLARE TAPE6=OUTPUT IN THE      C
!!!                  PROGRAM CARD).                                      C
!!!     IENT...      ENTRY PARAMETER FOR SUBROUTINE HISTO.               C
!!!                  FOR IENT=1, SUBROUTINE HISTO WILL BE INITIALIZED    C
!!!                  FOR PARAMETER NO. J.                                C
!!!                  FOR IENT=2, THE HISTOGRAM WILL BE UPDATED FOR       C
!!!                  SAMPLE VALUE X OF PARAMETER NO. J.                  C
!!!                  FOR IENT=3, THE HISTOGRAM WILL BE PRINTED FOR       C
!!!                  PARAMETER NO. J AND AGAIN INITIALIZED FOR PARAMETER C
!!!                  NO. J.                                              C
!!!     NAM...       ARRAY OF ALPHANUMERIC PARAMETER NAMES TO BE PRINTED C
!!!                  WITH FORMAT A10.                                    C
!!!     XMIN...      ARRAY OF LOWER BORDERS FOR HISTOGRAMS.              C
!!!     XMAX...      ARRAY OF UPPER BORDERS FOR HISTOGRAMS.              C
!!!     RDX...       ARRAY OF INCREMENTS FOR HISTOGRAMS.                 C
!!!     J...         PARAMETER NO.                                       C
!!!     X...         SAMPLE VALUE OF PARAMETER NO. J.                    C
!!!     NEWP...      PARAMETER FOR BEGINNING A NEW PAGE.                 C
!!!                  NEWP=1...HISTOGRAM WILL BE PRINTED ON A NEW PAGE    C
!!!                  FOR ENTRY PARAMETER IENT=3.                         C
!!!                  X HAS TO BE DEFINED ONLY FOR ENTRY PARAMETER IENT=2.C
!!!                                                                      C
!!!     PROGRAM CREATION...    820420 BY G.WEBER,                        C
!!!                            INSTITUT FUER ERDMESSUNG,                 C
!!!                            UNIVERSITAET HANNOVER,                    C
!!!                            NIENBURGER STR.6,                         C
!!!                            D-3000 HANNOVER 1,                        C
!!!                            FEDERAL REPUBLIC OF GERMANY.              C
!!!     PROGRAM MODIFICATION...830725 BY H.-G.WENZEL                     C
!!!     PROGRAM MODIFICATION...830909 BY H.-G.WENZEL                     C
!!!     PROGRAM MODIFICATION...890302 BY H. DENKER (FOR IBM)             C
!!!                                                                      C
!!!**********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER NAM*41, IZEI*1
      DIMENSION NAM(30),XMIN(30),XMAX(30),RDX(30),KMAX(30)
      DIMENSION RMX(30),RMS(30),RXMIN(30),RXMAX(30)
      DIMENSION N(52,30),IZEI(40)
      SAVE
      COMMON/HIS/ N
!!!     LEVEL2,/HIS/
      IF(IENT.NE.2) GOTO 100
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!     IENT=2. UPDATING OF ARRAYS N,RMX,RMS FOR PARAMETER NO. J         C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      K=(X-XMIN(J))/RDX(J)+2.
      IF(K.LT.1) K=1
      IF(K.GT.KMAX(J)) K=KMAX(J)
      N(K,J)=N(K,J)+1
      RMX(J)=RMX(J)+X
      RMS(J)=RMS(J)+X*X
      IF(X.LT.RXMIN(J)) RXMIN(J)=X
      IF(X.GT.RXMAX(J)) RXMAX(J)=X
      RETURN
  100 IF(IENT.NE.1) GOTO 300
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!     IENT=1. INITIALIZATION OF SUBROUTINE HISTO FOR PARAMETER NO. J.  C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 110 K=1,52
  110 N(K,J)=0
      RMX(J)=0.
      RMS(J)=0.
      RXMIN(J)= 999999.9
      RXMAX(J)=-999999.9
      KMAX(J)=(XMAX(J)-XMIN(J))/RDX(J)+2.
      IF(KMAX(J).GT.52) KMAX(J)=52
      RETURN
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!     IENT=3. PRINT HISTOGRAM FOR PARAMETER NO. J.                     C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  300 NSU=0
      NMA=0
      DO 310 K=1,KMAX(J)
      IF(N(K,J).GT.NMA) NMA=N(K,J)
  310 NSU=NSU+N(K,J)
      IF(NSU.EQ.0) RETURN
      SCALE=40./REAL(NMA)
      IF(NEWP.EQ.1) WRITE(IUN6,7006)
      WRITE(IUN6,7000)
      RMX(J)=RMX(J)/REAL(NSU)
      RRMS=SQRT(RMS(J)/REAL(NSU))
      RMS(J)=RMS(J)-RMX(J)*RMX(J)*REAL(NSU)
      RMS(J)=SQRT(ABS(RMS(J))/REAL(NSU))
      WRITE(IUN6,7003) NAM(J),NSU,RMX(J),RMS(J),RRMS,RXMIN(J),RXMAX(J)
      DO 320 K=1,KMAX(J)
      DO 330 L=1,40
  330 IZEI(L)=' '
      IZEI(1)='.'
      IZEIE=REAL(N(K,J))*SCALE
      IF(IZEIE.LE.40) GOTO 340
      IZEI(40)='+'
      IZEIE=40
  340 CONTINUE
      IF(IZEIE.LT.1) GOTO 360
      DO 350 L=1,IZEIE
  350 IZEI(L)='*'
  360 CONTINUE
      XMI=REAL(K-2)*RDX(J)+XMIN(J)
      XMA=REAL(K-1)*RDX(J)+XMIN(J)
      PROZ=REAL(N(K,J))/REAL(NSU)
      IF(K.NE.1.AND.K.NE.KMAX(J)) WRITE(IUN6,7001) K-1,XMI,XMA,PROZ,&
       N(K,J),(IZEI(L),L=1,40)
      IF(K.EQ.1) WRITE(IUN6,7002) K-1,XMA,PROZ,N(K,J),(IZEI(L),L=1,40)
      IF(K.EQ.KMAX(J)) WRITE(IUN6,7004) K-1,XMI,PROZ,N(K,J),&
       (IZEI(L),L=1,40)
  320 CONTINUE
  370 WRITE(IUN6,7005)
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!     RE-INITIALIZATION FOR PARAMETER NO.J.                            C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 380 K=1,52
  380 N(K,J)=0
      RMX(J)=0.
      RMS(J)=0.
      RXMIN(J)=999999.9
      RXMAX(J)=-999999.9
      RETURN
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!     FORMAT STATEMENTS.                                               C
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7000 FORMAT(/' SUBROUTINE HISTO, VERSION 03/02/89.'/&
              ' ===================================')
 7001 FORMAT(1X,I3,F10.3,' ... ',F10.3,F10.3,I10,3X,40A1)
 7002 FORMAT(1X,I3,' -INFINITY ... ',F10.3,F10.3,I10,3X,60A1)
 7004 FORMAT(1X,I3,F10.3,' ...   INFINITY',F10.3,I10,3X,40A1)
 7003 FORMAT(&
      ' HISTOGRAM FOR PARAMETER     ',10X,A10/&
      ' NUMBER OF SAMPLE VALUES     ',10X,I10/&
      ' MEAN OF SAMPLE VALUES       ',10X,F10.3/&
      ' STD. DEV. OF SAMPLE VALUES  ',10X,F10.3/&
      ' RMS OF SAMPLE VALUES        ',10X,F10.3/&
      ' MINIMUM SAMPLE VALUE        ',10X,F10.3/&
      ' MAXIMUM SAMPLE VALUE        ',10X,F10.3/)
 7005 FORMAT(/' END OF SUBROUTINE HISTO')
 7006 FORMAT(1H1)
      END SUBROUTINE

      endmodule


!!! program main

subroutine main1
    use intptFunc
    implicit none
    real*8, allocatable :: tElapsedT, VecFlat(:), VecFlon(:), VecU(:)
    real*8 :: VecFlatTmp(3), VecFlonTmp(3)

    allocate(VecFlat(3), VecFlon(3), VecU(3))

    data VecFlatTmp /0.0, 89.0,  -89.0/
    data VecFlonTmp /0.0, 179.0, 359.0/
    VecFlat = VecFlatTmp
    VecFlon = VecFlonTmp

    call INTPT(VecFlat, VecFlon, VecU, 3)

    print *, 'the main function end.'
end subroutine


