      SUBROUTINE NFERT(APMU,IRC,JFT,JRT)
!     APEX0806
!     THIS SUBPROGRAM APPLIES N AND P FERTILIZER AT SPECIFIED DATES,
!     RATES, AND DEPTH OR AUTOMATICALLY.
      USE PARM
      JRT=0
!     IF(ANA(IRO(ISA),ISA)>=FNMX(IRO(ISA),ISA).AND.IRC/=1)RETURN
      X1=0.
      I=LID(1,ISA)
      ZFT=TLD(JFT)
      SELECT CASE(IRC)
	      CASE(1,2)
	          I1=IABS(IAPL(ISA))
              KF=IDFT(IRC,I1)
	      CASE(3)
	          KF=IDMU(IHD,IOW)
	      CASE(4,5)
	          KF=IDFT(IRC,ISA)
	      CASE(6)
	          KF=LFT(IRO(ISA),KF1(ISA),ISA)
              X1=WFA(IRO(ISA),KF1(ISA),ISA)
              GO TO 5
          CASE(7)
              KF=IDFT(3,ISA)
              X1=APMU
	          X3=0.
	          GO TO 12              
      END SELECT
      IF(IRC/=3.AND.IRC/=4.AND.IRC/=7.AND.MNUL/=0)THEN
          XX=FP(KF)+FPO(KF)
          XX1=XMAP(ISA)/XX
          XX2=SAMA(ISA)/XX
          IF(XX2>=XX1)THEN
              JRT=1
              RETURN 
          END IF
          APMU=MIN(APMU,XX1-XX2)
      END IF
      X1=APMU
    5 IF(FOC(KF)<1.E-10)THEN
          DO J=1,NBSL(ISA)
              I=LID(J,ISA)
              IF(ZFT<Z(I,ISA))EXIT
          END DO
          IF(X1<1.E-10)THEN
              X3=MAX(UNA(JJK,ISA)-TNOR(ISA),0.)
              X1=X3/FN(KF)
              IF(X1>0.)GO TO 1
              IF(IRC==6)KF1(ISA)=KF1(ISA)+1
              RETURN
          END IF
      END IF
      X3=X1*FN(KF)
      GO TO 12
    1 X8=FNMX(IRO(ISA),ISA)-ANA(IRO(ISA),ISA)
      IF(X8<1.E-10)THEN
          JRT=1
          RETURN 
      END IF
      IF(X3<=X8)GO TO 12
      X3=X8
      X1=X3/FN(KF)
   12 X2=X1*FP(KF)
      X4=X3*FNMA(KF)
      X5=X1*FNO(KF)
      X6=X1*FPO(KF)
      X7=X3-X4
      X8=X1*FOC(KF)
      X11=X1*FSLT(KF)
      SMM(99,MO,ISA)=SMM(99,MO,ISA)+X8
      VAR(99,ISA)=X8
      IF(X8<.1)THEN
          WPML(I,ISA)=WPML(I,ISA)+X2
          WNMN(I,ISA)=WNMN(I,ISA)+X7
      ELSE
          RLN=.175*X8/(X3+X5)
          X10=.85-.018*RLN
          IF(X10<.01)THEN
              X10=.01
          ELSE
              IF(X10>.7)X10=.7
          END IF   
          XX=X8*X10
          WLMC(I,ISA)=WLMC(I,ISA)+XX
          YY=X1*X10
          WLM(I,ISA)=WLM(I,ISA)+YY
          !ZZ=X5*X10
          !WLMN(I,ISA)=WLMN(I,ISA)+ZZ
          !WLSN(I,ISA)=WLSN(I,ISA)+X5-ZZ
          XZ=X8-XX
          WLSC(I,ISA)=WLSC(I,ISA)+XZ
          WLSLC(I,ISA)=WLSLC(I,ISA)+XZ*.175
          WLSLNC(I,ISA)=WLSC(I,ISA)-WLSLC(I,ISA)
          YZ=X1-YY
          WLS(I,ISA)=WLS(I,ISA)+YZ
          WLSL(I,ISA)=WLSL(I,ISA)+YZ*.175
          WPMU(I,ISA)=WPMU(I,ISA)+X2
          WNMU(I,ISA)=WNMU(I,ISA)+X7
          Z1=.001*X1
          RSDM(I,ISA)=RSDM(I,ISA)+Z1
          WCOU(I,ISA)=WCOU(I,ISA)+X8
          WNOU(I,ISA)=WNOU(I,ISA)+X5
          WPOU(I,ISA)=WPOU(I,ISA)+X6
          SMM(78,MO,ISA)=SMM(78,MO,ISA)+Z1
          SMM(53,MO,ISA)=SMM(53,MO,ISA)+X5
          SMM(56,MO,ISA)=SMM(56,MO,ISA)+X6
      END IF
      ANA(IRO(ISA),ISA)=ANA(IRO(ISA),ISA)+X3
      WNMA(I,ISA)=WNMA(I,ISA)+X4
      SMM(54,MO,ISA)=SMM(54,MO,ISA)+X7
      SMM(55,MO,ISA)=SMM(55,MO,ISA)+X4
      SMM(57,MO,ISA)=SMM(57,MO,ISA)+X2
      IF(IRC>2)THEN
          TFNO=TFNO+X5*WSA(ISA)
          TFMN=TFMN+X7*WSA(ISA)
          TFMA=TFMA+X4*WSA(ISA)
      END IF
      VAR(53,ISA)=X5
      VAR(56,ISA)=X6
      VAR(54,ISA)=X7
      VAR(55,ISA)=X4
      VAR(57,ISA)=X2
      IF(ZFT<.01)THEN
          FSFN(ISA)=FSFN(ISA)+X5+X7+X4
          FSFP(ISA)=FSFP(ISA)+X6+X2
      END IF
      WSLT(I,ISA)=WSLT(I,ISA)+X11
      SMM(132,MO,ISA)=SMM(132,MO,ISA)+X11
      IF(IRC/=3)THEN
          FRTN(JJK,ISA)=FRTN(JJK,ISA)+X3+X5
          FRTP(JJK,ISA)=FRTP(JJK,ISA)+X2+X6
      END IF
      VAR(132,ISA)=X11
	  XX=X1*FCST(KF)
      COST(ISA)=COST(ISA)+XX
      IF(IRC==4)THEN
          Y1=COTL(JFT)
          Y2=Y1-COOP(JFT)
          COST(ISA)=COST(ISA)+Y1
          CSFX=CSFX+Y2
      END IF
      IF(KFL(31)>0)THEN
          WRITE(KW(31),34)ISA,NBSA(ISA),IYR,MO,KDA,FTNM(KF),KDC(JJK),&
          KDF(KF),IHC(JFT),NBE(JFT),NBT(JFT),XX,XX,X1
          IF(IRC==4)WRITE(KW(31),50)ISA,NBSA(ISA),IYR,MO,KDA,TIL(JFT),&
          KDC(JJK),IHC(JFT),NBE(JFT),NBT(JFT),Y1,Y2,FULU(JFT)
      END IF
      IF(NOP>0.OR.NBSA(ISA)==ISAP.AND.X3>0.)THEN
          IF(IRC>1)THEN
              WRITE(KW(1),11)ISA,NBSA(ISA),IYR,MO,KDA,IDON(ISA),IHD,FTNM(KF)&
              ,X1,ZFT,X3,X4,X5,X2,X6,XHSM(ISA)
          ELSE     
              WRITE(KW(1),25)ISA,NBSA(ISA),IYR,MO,KDA,IDON(ISA),NBSA(I1),&
              FTNM(KF),X1,ZFT,X3,X4,X5,X2,X6,XHSM(ISA)
          END IF   
      END IF
      IF(KFL(19)>0)WRITE(KW(19),22)ISA,NBSA(ISA),IYR,MO,KDA,IY,IDON&
      (ISA),IHD,FTNM(KF),X1,X3,X4,X5,X2,X6
      APMU=X1*WSA(ISA)
      IF(IRC==6)THEN
          KF1(ISA)=KF1(ISA)+1
          IF(X8<.1)GO TO 15
      END IF
      IF(IRC==4.OR.IRC==7)GO TO 15
      IF(NHRD(IOW)>0)THEN
          IF(IAPL(ISA)==0)RETURN
      END IF
      SMM(81,MO,ISA)=SMM(81,MO,ISA)+X1
      SAMA(ISA)=SAMA(ISA)+X2+X6
      TMAP=TMAP+APMU
	  OMAP(IDON(ISA))=OMAP(IDON(ISA))+APMU
!	  IF(IDON(ISA)==32)THEN
!	      WRITE(KW(1),101)ISA,NBSA(ISA),IY,MO,KDA,APMU,SMNU(32),&
!	      &XMAP(ISA),SAMA(ISA)
!	  END IF
! 101 FORMAT(1X,'!!!!!',2I8,3I4,4E16.6)	  
	  IF(IRC==3)RETURN
      GO TO 16
   15 FCMN(ISA)=FCMN(ISA)+X3+X5
      FCMP(ISA)=FCMP(ISA)+X2+X6
   16 NDFA(ISA)=0
      RETURN
   11 FORMAT(1X,2I8,1X,I4,2I2,2X,'IDON=',I4,2X,'HRD#=',I3,2X,A8,2X,&
      'RATE=',F6.0,'kg/ha',1X,'DPTH=',F5.2,'M ELEM WT(kg/ha)--',1X,'MN='&
      ,F4.0,1X,'NH3=',F4.0,1X,'ON=',F4.0,1X,'MP=',F4.0,1X,'OP=',F4.0,1X,&
      'HUSC=',F5.2)
   22 FORMAT(1X,2I8,1X,I4,2I2,3I4,3X,A8,F8.0,6F8.2)
   25 FORMAT(1X,2I8,1X,I4,2I2,2X,'IDON=',I4,2X,'FA#=',I4,2X,A8,2X,&
      'RATE=',F6.0,'kg/ha',1X,'DPTH=',F5.2,'M ELEM WT(kg/ha)--',1X,'MN='&
      ,F4.0,1X,'NH3=',F4.0,1X,'ON=',F4.0,1X,'MP=',F4.0,1X,'OP=',F4.0,1X,&
      'HUSC=',F5.2)
   34 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,2X,4I4,F10.2,10X,2F10.2)
   50 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
      END