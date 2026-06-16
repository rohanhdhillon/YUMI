SUBROUTINE BESSEL(L,Z,BJ,DBJ,BN,DBN)                              
        use mod_constantes, only : dp
implicit real(kind=dp) (a-h,o-z)
real(kind=dp), intent(out) :: bj,dbj,bn,dbn
real(kind=dp), intent(in) :: z
integer, intent(in) :: l
      IF(Z.GT.0.d0) GO TO 20                                             
      LP = L+1                                                          
      LM = L-1                                                          
      CALL BESM3(L,Z,G)                                                 
      BN = G                                                            
      CALL BESM3(LM,Z,G)                                                
      GM = G                                                            
      CALL BESM3(LP,Z,G)                                                
      GP = G                                                            
      DBN = -(L*GM+LP*GP)/(L+LP)                                        
      BJ = 0.d0                                                         
      DBJ= 0.d0                                                         
      RETURN                                                            
   20 CONTINUE                                                          
      CALL BESJOT(L,Z,F,DF,R)                                           
      CALL BESSEN(L,Z,G,DG,RG)                                          
      BJ = F                                                            
      DBJ = DF                                                          
      BN = G                                                            
      DBN = DG                                                          
      IF(L.GT.0)  GOTO 10                                               
      RETURN                                                            
 10   DO 1 I=1,L                                                        
      BJ = BJ / R                                                       
      DBJ = DBJ / R                                                     
      BN = BN * RG                                                      
    1 DBN = DBN * RG                                                    
      RETURN                                                            
      END                                                               
SUBROUTINE BESJOT(L,X,F,DF,R)                                     
!                                                                       
!   THIS PROGRAM GENERATES THE STANDARD AND MODIFIED VERSIONS OF THE SPH
!   BESSEL-FUNCTIONS OF FIRST(BESJOT)- AND SECOND(BESSEN)-KIND RESPECTIV
!   FOR DEFINITIONS COMPARE: 'NBS-H1NDBOOK OF MATHEMATICAL FUNCTIONS',  
!   (ABRAMOWITZ+STEGUN,EDS./N.Y.:1964), SECTIONS 10.1.1 ON PAGE 437 FOR 
!   STANDARD VERSIONS AND SS.10.2.2 + 10.2.3 ON P.443 FOR THE MODIFIED O
!   L=INDEX(NATURAL NUMBERS INCLUDING ZERO), X=ARGUMENT(REAL,D.P.), F=OU
!                                                                       
!   THE SIGN OF THE ARGUMENT IS USED TO DETERMINE THE VERSIONS:         
!   THE OUTCOMES F(=FIRST-KIND-FUNCTIONS) AND G(=SECOND-KIND-F.) MUST BE
!   DIVIDED (RESP. MULTIPLIED) BY THE L-TH POWER OF THE REDUCTION FACTOR
!   TO GET THE MODIFIED VERSIONS, USE ARGUMENT WITH NEGATIVE SIGN }     
!                                                                       
!   BY FORMULAS 10.1.31 ON PAGE 439 LOC.CIT. AND 10.2.7 ON P.443 IBID., 
!   SOLUTIONS HAVE BEEN TESTED TO BE CORRECT TO TWELVE PLACES AT LEAST I
!   RANGE COMBINING X=1...441 AND L=0...340 .                           
!                                                                       
!   BESJOT IS DIVIDED INTO THREE PARTS, CORRESPONDING TO WETHER X > L, O
!   WHILE X < L, BEEING 0.5*X*X < 2*L OR 0.5*X*X > 2*L RESPECTIVELY .   
!                                                                       
!                                                                       
! MODIF. POUR X PLUS GRAND QUE L DANS BESJOT  : R=1                     
!        QQ SOIT X DANS BESSEN : R=1                                    
!   POUR EVITER LES OVERFLOWS OU UNDERFLOWS DANS LE PROG. APPELE POUR 50
!      VERSION JAN. 77                                                  
!                                                                       
        use mod_constantes, only : dp
implicit real(kind=dp) (a-h,o-z)
real(kind=dp), intent(out) :: f,df,r
real(kind=dp), intent(in) :: x
integer, intent(in) :: l
!                                                                       
      F = 0.D0                                                          
      R = 1.D0                                                          
      IF(X)  51,50,52                                                   
 50   IF(L.EQ.0)  F=1.D0                                                
      DF = 0.D0                                                         
      RETURN                                                            
!                                                                       
 51   SINIX = DSINH(-X)                                                 
      COSIX = DCOSH(-X)                                                 
       W = -1.D0                                                        
      GOTO 53                                                           
 52   PI = 6.283185307179586D0                                          
      XR = DMOD(X,PI)                                                   
      SINIX = DSIN(XR)                                                  
      COSIX = DCOS(XR)                                                  
       W = +1.D0                                                        
 53    Z = 1.D0 / DABS(X)                                               
       A = DFLOAT(L)                                                    
       R = A * Z                                                        
      IF(DABS(X)-A) 2,2,1                                               
   2  IF(0.5D0*X*X-2.D0*A) 4,4,3                                        
!                                                                       
!   FOR THE FOLLOWING VERSION SEE PAGE 439, SECTION 10.1.19 LOC.CIT., AN
!   SECTION 10.1.11 ON PAGE 438 IBIDEM}                                 
!   FOR THE MODIFIED CASE LOOK UP SS.10.2.18 AND 10.2.13 ON PS.444 AND 4
!   THIS VERSION IS USED, IF X > L .                                    
!                                                                       
  1   F0 = SINIX * Z                                                    
      G0 =-COSIX * Z                                                    
        R=1.D0                                                          
      IF(L) 11,11,12                                                    
 11    F = F0                                                           
      DF =-G0 - F0*Z                                                    
      RETURN                                                            
 12   IF(L-1) 13,13,14                                                  
 13    F = W * (F0-COSIX) * Z                                           
      DF = F0   - 2.0D0*F*Z                                             
      RETURN                                                            
 14   F1 = W * (F0-COSIX) * Z                                           
      IF(L.EQ.2)  GOTO 15                                               
       J = L-2                                                          
      DO 10 I=1,J                                                       
      F2 = W * (F1*DFLOAT(2*I+1)*Z - F0  )                              
      F0 = F1                                                           
   10 F1 = F2                                                           
   15  F = W * (F1*DFLOAT(2*L-1)*Z - F0  )                              
      DF = F1   - dfloat(L+1)*F*Z                                             
      RETURN                                                            
!                                                                       
!   FOR THE FOLLOWING VERSION SEE PAGE 453, EXAMPLE 2 LOC.CIT.}         
!   THIS VERSION IS USED, IF X < L AND 0.5*X*X > 2*L .                  
!                                                                       
  3    N = A + 25.D0 + DSQRT(A)                                         
      B0 = 0.D0                                                         
      B1 = 1.D0                                                         
      DO 20 J=1,N                                                       
      B2 = W * (B1*DFLOAT(2*(N-J)+3)*Z - B0/R) / R                      
      B0 = B1                                                           
      IF(N-L-J) 22,21,22                                                
 21    F = B2                                                           
      GOTO 20                                                           
 22   IF(N-L+1-J)  20,23,20                                             
 23   DF = B2                                                           
   20 B1 = B2                                                           
      DF = W**(L-1) * (DF/B1) * SINIX*Z                                 
       F = W**L * (F/B1) * SINIX*Z                                      
      DF = DF*R - dfloat(L+1)*F*Z                                             
      RETURN                                                            
!                                                                       
!   FOR THE FOLLOWING VERSION SEE P1GE 437, FORMULA 10.1.2  LOC.CIT.}   
!   FOR THE MODIFIED CASE FORMULA 10.2.5 ON PAGE 443 IS VALID.          
!   THIS VERSION IS USED, IF X < L AND 0.5*X*X < 2*L .                  
!                                                                       
  4    Y = -W * 0.5D0 * X * X                                           
      S0 = DFLOAT(2*L-1)                                                
      S1 = DFLOAT(2*L+1)                                                
      P0 = 1.D0                                                         
      P1 = 1.D0                                                         
      C0 = 1.D0                                                         
      C1 = 1.D0                                                         
      DO 30 I=1,15                                                      
      S0 = S0 + 2.D0                                                    
      S1 = S1 + 2.D0                                                    
      P0 = Y*P0/(S0*DFLOAT(I))                                          
      P1 = Y*P1/(S1*DFLOAT(I))                                          
      C0 = C0 + P0                                                      
   30 C1 = C1 + P1                                                      
       Q = 1.D0                                                         
      IF(L.EQ.1)  GOTO 32                                               
       J = L - 1                                                        
      DO 31 I=1,J                                                       
   31  Q = Q * A/DFLOAT(2*I+1)                                          
       F = Q * A/DFLOAT(2*L+1) * C1                                     
      DF = Q*C0*R - DFLOAT(L+1)*F*Z                                     
      RETURN                                                            
  32   F = C1 / 3.D0                                                    
      DF = C0*R - 2.D0*F*Z                                              
      RETURN                                                            
!                                                                       
!                                                                       
           ENTRY BESSEN(L,X,G,DG,R)                                     
!                                                                       
!   SPHERICAL BESSEL-(AND MODIFIED BESSEL-) FUNCTIONS OF THE SECOND KIND
!   THIS VERSION IS VALID FOR ALL INDICES AND ARGUMENTS.                
!                                                                       
       G = 0.D0                                                         
       A = DFLOAT(L)                                                    
       R = 1.D0                                                         
      IF(X) 61,60,62                                                    
 60   WRITE(6,300)                                                      
!300   FORMAT(1H0,'*******  ARGUMENT OF SPHERICAL BESSEL-FUNCTION OF SECOND KIND SHOULD NOT BE ZERO }')                                    
 300   format('zero')
      RETURN                                                            
 61   SINIX = DSINH(-X)                                                 
      COSIX = DCOSH(-X)                                                 
       W = -1.D0                                                        
      GOTO 63                                                           
 62   PI = 6.283185307179586D0                                          
      XR = DMOD(X,PI)                                                   
      SINIX = DSIN(XR)                                                  
      COSIX = DCOS(XR)                                                  
       W = +1.D0                                                        
 63    Z = 1.D0 / DABS(X)                                               
      G0 =-W * COSIX * Z                                                
      F0 = W * SINIX * Z                                                
      IF(L) 41,41,42                                                    
 41    G = G0                                                           
      DG = W*F0 - G0*Z                                                  
      RETURN                                                            
 42        R=1.D0                                                       
       IF(L-1) 43,43,44                                                 
 43    G = W * (G0-SINIX)*Z                                             
      DG = G0   - 2.0D0*G*Z                                             
      RETURN                                                            
 44   G1 = W * (G0-SINIX)*Z                                             
      IF(L.EQ.2)  GOTO 45                                               
       J = L-2                                                          
      DO 40 I=1,J                                                       
      G2 = W * (G1*DFLOAT(2*I+1)*Z - G0  )                              
      G0 = G1                                                           
   40 G1 = G2                                                           
  45   G = W * (G1*DFLOAT(2*L-1)*Z - G0  )                              
      DG = G1   - DFLOAT(L+1)*G*Z                                       
      RETURN                                                            
      END                                                               
SUBROUTINE BESM3(L,Z,G)                                           
!                                                                       
!   CESOUS PROGRAMME DONNE LES FONCTIONS DE BESSEL MODIFIEES DE3EME ESPE
!    MULTIPLIEES PAR EXP(Z) POUR EVITER LES UNDER FLOWS AUX GRANDES VALE
!     ABS(Z)   , VERSION JAN. 77                                        
!                                                                       
        use mod_constantes, only : dp
implicit real(kind=dp) (a-h,o-z)
real(kind=dp), intent(out) :: g
real(kind=dp), intent(in) :: z
integer, intent(in) :: l
real(kind=dp) :: fct(1000)
      DATA PI /3.141592653589793D0/                                     
      FCT (1) = 0.D0                                                    
      DO  I = 1 , 999                                                  
      AI = dfloat(I)                                                            
      FCT (I+1) = FCT(I) + DLOG(AI)                                     
      enddo
      G = 0.d0                                                            
      IF(L.LT.0) RETURN                                                 
      ZM = DABS(Z)                                                      
      D=0.5D0/ZM                                                        
      S= D                                                              
      IF(L.EQ.0) GO TO 100                                              
      DL=DLOG(D)                                                        
      DL1=0.D0                                                          
      L1 = L+1                                                          
      S = 0.d0                                                            
      DO 1 I = 1 , L1                                                   
      K = I-1                                                           
      DL1 = DL + DL1                                                    
!     write(*,*)'L,K',L,K
      ARG=      FCT(L+K+1)-FCT(K+1)-FCT(L-K+1) +DL1                     
      IF (ARG.LT.-180.d0) GO TO 100                                       
      S= S+DEXP(ARG)                                                    
    1 CONTINUE                                                          
  100 G = PI*S                                                          
      RETURN                                                            
      END                                                               

