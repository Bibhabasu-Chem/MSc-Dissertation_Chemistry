C******************TRYRb.FOR*******************************
C*********CALCULATION OF 'Q' FROM S(k)***************
       SUBROUTINE SKMET (AK,SK)
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
       PI=3.142857143
       ETA=(PI*RHO*SIG**3)/6.0
       GAMAB=ETA*(1.0+2.0*ETA)**2/(2.0*(1.0-ETA)**4)
       BETAB=-6.0*ETA*(1.0+ETA*0.5)**2/(1.0-ETA)**4
       ALPHAB=(1.0+2.0*ETA)**2/(1.0-ETA)**4
       ASI=AK*SIG
       RCK=(-24.0*ETA/(ASI)**6)*(ALPHAB*(ASI)**3*(SIN(ASI)-ASI*
     & COS(ASI))+BETAB*(ASI)**2*(2.0*ASI*SIN(ASI)-((ASI)**2-2.0)*
     & COS(ASI)-2.0)+GAMAB*((4.0*(ASI)**3-24.0*ASI)*SIN(ASI)-
     & ((ASI)**4-12.0*(ASI)**2+24.0)*COS(ASI)+24.0)-(EPS*(ASI)**3/
     & TEMP)*(SIN(ALAM*ASI)-ALAM*ASI*COS(ALAM*ASI)+ASI*COS(ASI)-
     & SIN(ASI)))
       SK=1.0/(1.0-RCK)
C       SKP=SK*1.5
       RETURN
       END


       SUBROUTINE GRMET (R,GR)
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
       COMMON/SFS/S11(200)
       A(SI)=AK*SIN(AK*SI)
       GRI(SIJ,SG)=(SIJ-1.0)*A(SG)*EXP(-4.5*AK*AK/(13.0)**2)
C       GRI(SIJ,SG)=(SIJ-1.0)*A(SG)
       SG1=0.0
       PI=3.142857143
       H=0.1/3.0
       GD1=1.0/(2.0*PI*PI)
       DO 18 J=3,291
       AK=.1*FLOAT(J)
       IF (J.EQ.3.OR.J.EQ. 291)CZ=1.0
       AG1=GRI(S11(J),R)
       SG1=SG1+AG1*CZ
       IF (J.EQ.3)CZ=2.0
18     CZ=6.0-CZ
       SG1=SG1*H
       GS11=1.0+GD1/(RHO*R)*SG1
       GR=GS11
       RETURN
       END


       SUBROUTINE ZIE
       COMMON/SFS/S11(500)
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
       AM=65.38
       PI=3.142857143
       ETA=(PI*RHO*SIG**3)/6.0
       CONST=(8.0/3.0)*SQRT(PI*AM*1.38/0.6023)
C       CALL GRMET (2.5,GRSIG)
C       WRITE(*,*) 'g(sig)=  ', GRSIG
       GRSIG=(1.0+ETA/2.0)/(1.0-ETA)**2
C       WRITE(*,*)GRSIG
       ZIH=CONST*SIG**2*SQRT(TEMP)*GRSIG*RHO
       WRITE(6,95)ZIH
95     FORMAT(2X,'ZIH=',E11.4)
       RETURN
       END


       SUBROUTINE ZIES
       COMMON/SFS/S11(500)
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
       A(SI)=AK*SIN(AK*SI)
       VIJ(AII,SIGZ)=(AII*AK*SIGZ*COS(AII*AK*SIGZ)-
     & SIN(AII*AK*SIGZ)-AK*SIGZ*COS(AK*SIGZ)+SIN(AK*SIGZ))
C       GRI(SIJ)=(SIJ-1.0)*EXP(-1.5*AK*AK/(JZ5*.1)**2)
       GRI(SIJ)=(SIJ-1.0)
       FZI(SL,AII,SGI)=GRI(SL)*VIJ(AII,SGI)
       H=0.1/3.0
       PI=3.142857143
       CON=-(EPS*1.38)/3.0*SQRT(1.0/(PI*0.6023*1.38))
       AM=65.38
       TAM=SQRT(TEMP*AM)
       DO 180 JZ5=61,61
       S1=0.0
       DO 18 J=1,JZ5
       AK=0.1*FLOAT(J)
       IF (J.EQ.1 .OR. J.EQ.JZ5)CZ=1.0
       ASI1=FZI(S11(J),ALAM,SIG)
       S1=S1+ASI1*CZ
       IF (J.EQ.1)CZ=2.0
       CZ=6.0-CZ
18     CONTINUE
       S1=S1*H
       ZIS1=CON*TAM*S1
       ZIS1=ZIS1/TEMP
C       WRITE(6,97)JZ5,ZIS1,ASI1
       WRITE(6,97)ZIS1
180    CONTINUE
C97     FORMAT(2X,'JZ5=',I5,3X,'ZIS1=',E11.4, 6X, E11.4)
97     FORMAT(2X,'ZIS1=',E11.4)
       RETURN
       END


       SUBROUTINE ZIESH
       COMMON/SFS/S11(500)
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
       VIJ(AII,SIGZ)=(AII*AK*SIG*COS(AII*AK*SIG)-SIN(AII*AK*SIG)
     & -AK*SIG*COS(AK*SIG)+SIN(AK*SIG))
       GRI(SG)=(AK*SG*COS(AK*SG)-SIN(AK*SG))/(AK**3)
       FZI(SG,AII,SIGI)=GRI(SG)*VIJ(AII,SIGI)
       PI=3.142857143
       AM=65.38
       H=0.1/3.0
       SG=2.5
       SH1=0.0
       ETA=(PI*RHO*SIG**3)/6.0
       GRSIG=(1.0+ETA/2.0)/(1.0-ETA)**2
C       CALL GRMET(4.306,GRSIG)
C       WRITE(*,*)GRSIG
       CON=(-1.38*4.0/3.0)*EPS*RHO*GRSIG*SQRT(AM*PI/
     & (1.38*0.6023*TEMP))
       DO 18 J=1,61
       AK=0.1*FLOAT(J)
       IF (J. EQ. 1. OR. J. EQ. 61)CZ=1.0
       ASHI=FZI(SG,ALAM,SIG)
       SH1=SH1+ASHI*CZ
       IF (J. EQ. 1)CZ=2.0
       CZ=6.0-CZ
18     CONTINUE
       SH1=SH1*H
       ZISH=CON*SH1
       WRITE(6,99)ZISH
99     FORMAT(2X,'ZISH=',E11.4)
       RETURN
       END


       SUBROUTINE DIF
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
       ZIH=8.538
       ZIS1=1.369
       ZISH=1.753
       D=(1.38*TEMP/(ZIH+ZIS1+ZISH))*0.1
       WRITE(6,100)D
100    FORMAT(2X,'D=',E11.4)
       RETURN
       END


       SUBROUTINE DIFT
       COMMON/SFS/S11(500)
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
       PI=3.142857143
       ZIH=8.538
       ZIS1=1.369
       ZISH=1.753
       ZIE=ZIH+ZIS1+ZISH
       ETA=PI*RHO*SIG**3/6.0
       ALPHA=(1.0+2.0*ETA)**2/(1.0-ETA)**4
       BETA=-6.0*ETA*(1.0+ETA*0.5)**2/(1.0-ETA)**4
       GAMA=(ETA/2.0)*(1.0+2.0*ETA)**2/(1.0-ETA)**4
       AM=85.467
       BK=1.38
       H=.1/3.
       GSIG=(1.0+ETA/2.0)/(1.0-ETA)**2
C       CALL GRMET(2.5,GSIG)
C       WRITE(*,*)GSIG
C       ALP1=((1.0+ETA+ETA**2)/(1.0-ETA)**3)/((TEMP*(1.0+4.0*ETA**2
C     & +4.0*ETA)/(1.0-ETA)**4)-(8.0*EPS*ETA*(ALAM**3-1)))
C       WRITE(*,*)ALP1
       ALP1=0.000339
       CON1=(2.0/3.0)*ALP1*(GSIG-1.0)
       CON2=ALP1/(6.0*PI*PI*RHO)
       CON3=1.0/(2.0*PI*PI*RHO*SIG)
       CON4=(RHO/(12.0*PI*PI))*SQRT(PI*AM/(0.6023*BK*TEMP))
       CON5=(RHO/3.0)*SQRT(AM/(PI*0.6023*BK*TEMP))
       AIN1=0.0
       AIN2=0.0
       AIN3=0.0
       AIN4=0.0
       DO 12 J=1,61
       AK=0.1*FLOAT(J)
       ASG=AK*SIG
       AP=ASG**3*(SIN(ASG)-ASG*COS(ASG))
       AQ=ASG**2*(2.0*ASG*SIN(ASG)-(ASG**2-2.0)*COS(ASG)-2.0)
       AR=((4.0*ASG**3-24.0*ASG)*SIN(ASG)-(ASG**4-12.0*
     & ASG**2+24.0)*COS(ASG)+24.0)
       AS=ASG**3*(SIN(ALAM*ASG)-ALAM*ASG*COS(ALAM*ASG)+
     & ASG*COS(ASG)-SIN(ASG))
       ZK=ALP1*S11(J)*((1.0-S11(J))*EXP(-1.5*AK*AK/(6.1)**2)+
     & (24.0*ETA*S11(J)/(ASG**6))*(AP*4.0*ALPHA*ETA*(2.0+ETA)
     & /((2.0*ETA+1.0)*(1.0-ETA))+AQ*BETA*(ETA**2+9.0*ETA+2.0)
     & /((ETA+2.0)*(1.0-ETA))+AR*GAMA*(2.0*ETA**2+9.0*ETA+1.0)
     & /((1.0+2.0*ETA)*(1.0-ETA))-(EPS*AS)/(ALP1*TEMP**2)))
       GK=(S11(J)-1.0)*EXP(-1.5*AK*AK/(6.1)**2)/RHO
       FISK=(4.0*PI*EPS*BK/AK**3)*(ALAM*ASG*COS(ALAM*ASG)
     & -SIN(ALAM*ASG)-ASG*COS(ASG)+SIN(ASG))
       IF(J. EQ. 1. OR. J. EQ. 61)CZ=1.0
       AKS=((S11(J)-1.0)*EXP(-1.5*AK*AK/(6.1)**2))*AK**2*COS(ASG)
C       AKS=((S11(J)-1.0)*AK**2*COS(ASG))
       ZSK=AK*SIN(ASG)*ZK
       FZK=AK**3*FISK*(ZK/RHO+GK*ALP1)
       FSK=(ASG*COS(ASG)-SIN(ASG))*FISK
       AIN1=AIN1+AKS*CZ
       AIN2=AIN2+ZSK*CZ
       AIN3=AIN3+FZK*CZ
       AIN4=AIN4+FSK*CZ
       IF (J. EQ. 1)CZ=2.0
       CZ=6.0-CZ
12     CONTINUE
       AIN1=AIN1*H
       AIN2=AIN2*H
       AIN3=AIN3*H
       AIN4=AIN4*H
       DZIHT=ZIH/(2.0*TEMP)-ZIH*ALP1+8.0/3.0*RHO*SIG**2*
     & SQRT(PI*AM*BK*TEMP/0.6023)*(CON1+CON2*AIN1+CON3*AIN2)/100.0
       DZIST=(-ZIS1/(2.0*TEMP))-ALP1*ZIS1-(CON4*AIN3)/100.0
       DZISHT=(-ZISH/(2.0*TEMP))-ALP1*ZISH-CON5*AIN4*
     & (CON1+CON2*AIN1+CON3*AIN2)/100.0
       DZIT=DZIHT+DZIST+DZISHT
       DLNDT=(1.0/TEMP)-(1.0/ZIE)*DZIT
       AEQ=DLNDT*0.008314*TEMP*TEMP
       WRITE(6,101)DZIT,DLNDT,AEQ
C       WRITE(*,*)ALP1
101     FORMAT(2X,'DZIT=',F11.4,2X,'DLNDT=',F5.4,2X,'AEQ=',F11.7)
       RETURN
       END


C********************MAIN PROGRAM*************************
       COMMON/SFS/S11(200)
       COMMON/SRT/SIG,RHO,TEMP,ALAM,EPS
C       OPEN (UNIT=7,FILE='SKRb',STATUS='NEW')
C       OPEN (UNIT=6,FILE='DFT',STATUS='NEW')
C       OPEN (UNIT=6,FILE='ZIRb',STATUS='NEW')
C       OPEN (UNIT=6,FILE='ZIS',STATUS='NEW')
C       OPEN (UNIT=6,FILE='ZSH',STATUS='NEW')
C       OPEN (UNIT=6,FILE='GRRb',STATUS='NEW')
C      OPEN (UNIT=6,FILE='DIFEL',STATUS='NEW')
       PI=3.142857143   
       SIG=2.65
       RHO=0.04373 
       TEMP=443.
       ALAM=1.68
       EPS=200.0
       DO 105 J=1,200
       AK=0.1*FLOAT(J)
       CALL SKMET (AK,SK)
C       SKP=SK*1.5
        WRITE(*,110)AK
C      WRITE(7,110)AK,SKP
C160    FORMAT(2X,F5.4,',',2X,F10.5)
       S11(J)=SK
105    CONTINUE
       LL3=1000
       DO 25 JF=230,LL3
       R1=0.01*FLOAT(JF)
       WRITE(*,110) GR1
        CALL GRMET (R1,GR1)
25     CONTINUE
110    FORMAT(2X,F7.4,3X,F11.3)
       CALL ZIE
       CALL ZIES
       CALL ZIESH
       CALL DIF
       CALL DIFT
       STOP
       END

C  ZIH= 0.8538E+03 = 8.538E-10
C  ZIS1= 0.1369E+03 = 1.369E-10
C ZISH= 0.1753E+03 = 1.753E-10
C D= 0.3704E+01 = 3.707E-5
C  DZIT=     0.0077  DLNDT=.0025  AEQ=  2.0626776 KJ/MOLE
