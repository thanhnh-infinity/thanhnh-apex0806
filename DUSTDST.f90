      SUBROUTINE DUSTDST
!     APEX0806
!     THIS SUBPROGRAM SETTLES AND DISTRIBUTES PARTICULATE MATTER(DUST) 
!     FROM FEEDYARDS TO SURROUNDING SUBAREAS.
      USE PARM
      DIMENSION NX(100)
      DIMENSION DPY(100),WDRX(16)
      DATA A90/1.5708/,A180/3.1416/,A270/4.7124/,A360/6.2832/,WDRX/.1963&
     &,.589,.9817,1.374,1.767,2.16,2.553,2.945,3.338,3.731,4.123,4.516,&
     &4.909,5.301,5.694,6.087/
      IF(IDFH(ISA)==0)RETURN
!     IF(IFED(IDFH(ISA),IDON(ISA))/=ISA)RETURN
      TTH=TAN(THW)
      GRDX=.0005*GRDL(ISA)
      DO I=1,16
          IF(THW<WDRX(I))GO TO 4
      END DO
      I=1
    4 NWDR(I)=NWDR(I)+1
      CALL DUSTEM
      IF(DEMR<1.E-5)RETURN
!     Q=86400.*U10(IRF(ISA))*ZW*GRDL(ISA)
!     CW=PRMT(65)*U10(IRF(ISA))**PRMT(66)
!     CIN=DEMR/Q
!     DY=.1
!     X1=DEMR
!     TRT=DY/(3.6*U10(IRF(ISA)))
!     FY=1.-EXP(-PRMT(64)*TRT*10.)
!     SUM=0.
!     DO 25 I=1,40
!     SUM=SUM+DY
!     YDST(I)=SUM
!     DEP=Q*(CIN-CW)*FY
!     IF(DEP<0.)GO TO 26
!     DPY(I)=DEP
!     X1=X1-DEP
!     CIN=X1/Q
!  25 CONTINUE
!  26 NYD=I-1 
      SUM=0.
      NDX=0
      DO I=1,MSA
          DPY(I)=0.
          NX(I)=I
          IF(I==ISA)CYCLE
    !     COMPUTE BASE Y COORDINATE--X COORDINATE = WIND DIRECTION (THW)
          YB=(XCT(ISA)-XCT(I))*TTH+YCT(ISA)
    !     DETERMINE WHICH SUBAREAS ARE DOWN WIND FROM FEEDYARD
          IF(THW>A90.AND.THW<A270)GO TO 3
    !     CASE 1 OR 2 (N WIND)
          IF(YCT(I)>YB)CYCLE
          Y0=YCT(I)-YCT(ISA)
          X0=XCT(I)-XCT(ISA)
          DIST=SQRT(X0*X0+Y0*Y0)
          IF(THW<A90)GO TO 11
    !     CASE 2      
          IF(X0<0.)GO TO 12
          IF(Y0>0.)GO TO 13
          IF(Y0<X0*TTH)GO TO 14
    !     CASE 2.3
          ALF=ATAN(Y0/(X0+.001))
          BTA=THW-A270+ALF
          GO TO 10
    !     CASE 2.1
       12 ALF=ATAN(X0/(Y0+.001))
          BTA=THW-A270-ALF
          GO TO 9
    !     CASE 2.4
       13 ALF=ATAN(Y0/(X0+.001))
          BTA=A360-THW-ALF
          GO TO 9
    !     CASE 2.2
       14 ALF=ATAN(X0/(Y0+.001))
          BTA=A360-THW+ALF
          GO TO 10
    !     CASE 1
       11 IF(Y0>0.)GO TO 15
          IF(X0>0.)GO TO 16
          IF(Y0>X0/TTH)GO TO 17
    !     CASE 1.3
          ALF=ATAN(X0/(Y0+.001))
          BTA=THW-ALF
          GO TO 10
    !     CASE 1.1
       15 ALF=ATAN(Y0/(X0+.001))
          BTA=THW+ALF
          GO TO 9
    !     CASE 1.4
       16 ALF=ATAN(X0/(Y0+.001))
          BTA=A90-THW+ALF
          GO TO 9
    !     CASE 1.2
       17 ALF=ATAN(Y0/(X0+.001))
          BTA=A90-THW-ALF
          GO TO 10
        3 IF(YCT(I)<YB)CYCLE
    !     CASE 3 OR 4 (S WIND)
          Y0=YCT(I)-YCT(ISA)
          X0=XCT(I)-XCT(ISA)
          DIST=SQRT(X0*X0+Y0*Y0)
          IF(THW<A180)GO TO 18
    !     CASE 4
          IF(X0<0.)GO TO 22
          IF(Y0<0.)GO TO 23
          IF(Y0>X0/TTH)GO TO 24
    !     CASE 4.3
          ALF=ATAN(Y0/(X0+.001))
          BTA=A270-THW-ALF
          GO TO 10
    !     CASE 4.1
       22 ALF=ATAN(X0/(Y0+.001))
          BTA=A270-THW+ALF
          GO TO 9
    !     CASE 4.4
       23 ALF=ATAN(Y0/(X0+.001))
          BTA=THW-A180+ALF
          GO TO 9
    !     CASE 4.2
       24 ALF=ATAN(X0/(Y0+.001))
          BTA=THW-A180-ALF
          GO TO 10
    !     CASE 3
       18 IF(Y0<0.)GO TO 19
          IF(X0>0.)GO TO 20
          IF(Y0<X0*TTH)GO TO 21
    !     CASE 3.3
          ALF=ATAN(X0/(Y0+.001))
          BTA=A180-THW+ALF
          GO TO 10
    !     CASE 3.1
       19 ALF=ATAN(Y0/(X0+.001))
          BTA=A180-THW-ALF
          GO TO 9
    !     CASE 3.4
       20 ALF=ATAN(X0/(Y0+.001))
          BTA=THW-A90-ALF
          GO TO 9
    !     CASE 3.2
       21 ALF=ATAN(Y0/(X0+.001))
          BTA=THW-A90+ALF
          GO TO 10
    !     COMPUTE NEW X AND Y COORDINATES ON BASE (XB,YB)
        9 XT=ABS(DIST*COS(BTA))
          YT=ABS(DIST*SIN(BTA))
          GO TO 6
       10 XT=ABS(DIST*SIN(BTA))
          YT=ABS(DIST*COS(BTA))
        6 IF(XT<GRDX)THEN
              FX=1.
          ELSE
              X1=YT/XT
              IF(X1<.01)CYCLE
              X2=ATAN(X1)
              FX=(SIN(X2))**PRMT(67)
          END IF
          TRT=YT/(3.6*U10(IRF(ISA)))
          X1=PRMT(64)*TRT*10.
          IF(X1>10.)CYCLE
          FY=EXP(-X1)
          DPY(I)=FX*FY*WSA(I)
          SUM=SUM+DPY(I)
          NDX=NDX+1
    !     WRITE(KW(1),1)I,NBSA(I),IY,MO,KDA,THW,ALF,BTA,X0,Y0,XT,YT,DEMR,Q,
    !    1U10(IRF(ISA)),TRT,CIN,CW,FX,FY,DPY(I)
      END DO
      IF(SUM>0.)THEN
          B1=DEMR/SUM
      ELSE
          B1=1.
      END IF   
      TOT=0.
      DO I=1,MSA
          DPY(I)=DPY(I)*B1
          TOT=TOT+DPY(I)
          PM10(I)=PM10(I)+DPY(I)
          SMM(94,MO,I)=SMM(94,MO,I)+DPY(I)
      END DO
      CALL ASORT1(DPY,NX,MSA)
      WRITE(KW(20),32)IY,IYR,MO,KDA,THW,U10(IRF(ISA)),DEMR,TOT
      SUM=0.
      DO I=1,NDX
          X1=DPY(NX(I))/TOT
          SUM=SUM+X1
          IF(SUM>.999)CYCLE
    !     WRITE(KW(1),30)I,NX(I),NBSA(NX(I)),DPY(NX(I)),X1,SUM
          IF(KFL(20)>0)WRITE(KW(20),31)I,NX(I),NBSA(NX(I)),DPY(NX(I)),X1,&
         &SUM
      END DO
      !WRITE(KW(1),29)DEMR,TOT
      RETURN
!   1 FORMAT(1X,3I4,2I2,20E13.5)
!  29 FORMAT(1X,2E13.5)
!  30 FORMAT(1X,3I4,3E13.5)
   31 FORMAT(5X,I8,1X,2I8,F10.2,2F10.4)
   32 FORMAT(2X,'YR#=',I4,1X,'Y M D=',I4,2I2,1X,'W DIR=',F5.1,' RAD',2X,&
     &'W SPD=',F6.2,' m/s',2X,'DEMR=',F6.0,' KG',2X,'DDPR=',F6.0,' KG')
      END