      SUBROUTINE CGROW(JRT)
!     APEX0806
!     THIS SUBPROGRAM CALCUALTES LEAF AREA INDEX, HEAT UNITS, ROOT DEPTH
!     AND TEMPERATURE STRESS FOR THE CROP.
      USE PARM
!	  DIMENSION SLA0(MNC,MSA),XDLA0(MNC,MSA) 
      JRT=0
	  X1=DM(JJK,ISA)+1.E-10 
      CPR=UP1(JJK,ISA)/X1
      CNR=UN1(JJK,ISA)/X1
      AJWA=1.
      XPHU=PHU(JJK,IHU(JJK,ISA),ISA)
      TGX=TX-TBSC(JJK)
      IF(TGX>0.)HU(JJK,ISA)=HU(JJK,ISA)+TGX
      IF(JDA==JDHU.AND.IDC(JJK)/=NDC(7).AND.IDC(JJK)/=NDC(8).AND.IDC&
      (JJK)/=NDC(10))HU(JJK,ISA)=XPHU*PRMT(26)
      HUI(JJK,ISA)=HU(JJK,ISA)/XPHU
      IF(HU(JJK,ISA)>XPHU)THEN
          WCYD=MAX(WCY(JJK),WCYD-EO*.002)
          IF(IDC(JJK)==NDC(3).OR.IDC(JJK)==NDC(6))THEN
              JRT=2
              RETURN
          ELSE
              JRT=1
              RETURN
          END IF
      END IF
      IF(IDC(JJK)==NDC(8).OR.IDC(JJK)==NDC(10))THEN
          X1=HSM(ISA)/AHSM
          F1=X1/(X1+EXP(DLAP(1,JJK)-DLAP(2,JJK)*X1))
	      F=F1
          XLAI(JJK,ISA)=MAX(.1,DMLX(JJK)*HUI(JJK,ISA)/(HUI(JJK,ISA)+EXP(&
          DLAP(1,JJK)-DLAP(2,JJK)*HUI(JJK,ISA))))
      END IF
      F2=HUI(JJK,ISA)/(HUI(JJK,ISA)+EXP(DLAP(1,JJK)-DLAP(2,JJK)*&
      HUI(JJK,ISA)))
	  IF(IDC(JJK)/=NDC(8).AND.IDC(JJK)/=NDC(10))F=F2 
	  IF(IDC(JJK)==NDC(8).OR.IDC(JJK)==NDC(7).OR.IDC(JJK)==NDC(10))THEN
          F3=HUI(JJK,ISA)
	  ELSE
          F3=SQRT(F2+1.E-10) 
	  END IF
      FF=F-WLV(JJK,ISA)
      XX=FF*XLAI(JJK,ISA)
      X2=1.
      SLAX=0.
      X3=(SLAI(JJK,ISA)+.001)*CPHT(JJK,ISA)
      IF(IGO(ISA)==1)THEN
          SUM=X3
          TOT=STL(JJK,ISA)
          ADD=SLAI(JJK,ISA)
      ELSE
          SUM=0.
	      TOT=0.
	      ADD=0.
          DO I=1,IGO(ISA)
              K1=JE(I,ISA)
              IF(SLAI(K1,ISA)>SLAX)SLAX=SLAI(K1,ISA)
	          TOT=TOT+STL(K1,ISA)
              SUM=SUM+SLAI(K1,ISA)*CPHT(K1,ISA)
              ADD=ADD+SLAI(K1,ISA)
          END DO
          IF(SLAX>2.)X2=X3/SUM
	  END IF
      IF(XX>0.)THEN
          IF(IDC(JJK)/=NDC(7).AND.IDC(JJK)/=NDC(8).AND.IDC(JJK)/=NDC(10))THEN
              X1=X1*SQRT(REG(JJK,ISA))*SHRL
          ELSE
              X1=XX*X2*(1.+HR1)**3
          END IF
          SLAI(JJK,ISA)=MIN(XLAI(JJK,ISA),SLAI(JJK,ISA)+X1)
      END IF
      WLV(JJK,ISA)=F
      X1=TOPC(JJK)-TBSC(JJK)
      RTO=TGX/X1
      IF(RTO<2..AND.TGX>0.)THEN
          REG(JJK,ISA)=SIN(1.5707*RTO)
      ELSE
          REG(JJK,ISA)=0.
      END IF
      IF(SLAI(JJK,ISA)>.05)THEN
          CPHT(JJK,ISA)=MAX(CPHT(JJK,ISA),HMX(JJK)*F3)
          IF(HUI(JJK,ISA)<=XDLAI(JJK).OR.HRLT<=HLMN(ISA))THEN
              XDLA0(JJK,ISA)=1.-XDLAI(JJK)
              SLA0(JJK,ISA)=SLAI(JJK,ISA)
              GO TO 8
          END IF
          XX=MAX(1.E-10,(1.-HUI(JJK,ISA))/XDLA0(JJK,ISA))
          XX=LOG10(XX)
          IF(IDC(JJK)/=NDC(7).AND.IDC(JJK)/=NDC(8).AND.IDC(JJK)/=NDC(10))THEN
              RTO=RLAD(JJK)*XX
              IF(RTO<-10.)RTO=-10.
              SLAI(JJK,ISA)=MIN(SLAI(JJK,ISA),SLA0(JJK,ISA)*10.**RTO)
          END IF
          RTO=RBMD(JJK)*XX
          IF(RTO<-10.)RTO=-10.
          AJWA=10.**RTO
          GO TO 8
      END IF
      SLAI(JJK,ISA)=.05
    8 RD(JJK,ISA)=MIN(RDMX(JJK),Z(LID(NBSL(ISA),ISA),ISA),2.5*RDMX&
      (JJK)*HUI(JJK,ISA))
	  FGC(ISA)=ADD/(ADD+EXP(SCRP(23,1)-SCRP(23,2)*ADD)) 
	  FGSL(ISA)=TOT/(TOT+EXP(SCRP(24,1)-SCRP(24,2)*TOT)) 
      CLG(ISA)=BLG(3,JJK)*HUI(JJK,ISA)/(HUI(JJK,ISA)+EXP(BLG(1,JJK)-BLG&
      (2,JJK)*HUI(JJK,ISA)))
      SHRL=1.
      FHR=0.
      IF(HRLT+1.E-5<WDRM(ISA))THEN
          SHRL=0.
          FHR=1.-HRLT/WDRM(ISA)
      ELSE
          SRA(JJK,ISA)=SRA(JJK,ISA)+SRAD(IRF(ISA))
      END IF
      IF(IDC(JJK)==NDC(7))GO TO 14
      IF(TMN(IRF(ISA))>-1.)GO TO 11
      XX=ABS(TMN(IRF(ISA)))
      F=XX/(XX+EXP(FRST(1,JJK)-FRST(2,JJK)*XX))
      F=MAX(F,FHR)
      GO TO 12
   11 IF(SHRL>0.)GO TO 14
      F=FHR
   12 IF(IDC(JJK)/=NDC(8).AND.IDC(JJK)/=NDC(10))THEN
          IF(STL(JJK,ISA)>.01)THEN
              XX=F*STL(JJK,ISA)
              STL(JJK,ISA)=STL(JJK,ISA)-XX
              DM(JJK,ISA)=DM(JJK,ISA)-XX
              STD(JJK,ISA)=STD(JJK,ISA)+XX
              STDL(JJK,ISA)=STDL(JJK,ISA)+CLG(ISA)*XX
              XY=XX*CNR
              XZ=XX*CPR
              XUN=UN1(JJK,ISA)
              IF(XUN-XY<.01)XY=XUN-.01
              XUP=UP1(JJK,ISA)
              IF(XUP-XZ<.01)XZ=XUP-.01
              STDN(JJK,ISA)=STDN(JJK,ISA)+XY
              STDP(JJK,ISA)=STDP(JJK,ISA)+XZ
              UN1(JJK,ISA)=XUN-XY
              UP1(JJK,ISA)=XUP-XZ
          END IF
          IF(F*(1.-SNOF)>.9)THEN
              IF(IDC(JJK)/=NDC(3).AND.IDC(JJK)/=NDC(6))IGO(ISA)=IGO(ISA)-1
          END IF
      END IF
      SLAI(JJK,ISA)=SLAI(JJK,ISA)*(1.-F)
   14 IF(REG(JJK,ISA)>0.)RETURN
      SFCP(4,JJK,ISA)=SFCP(4,JJK,ISA)+1.
      SFMO(4,JJK,ISA)=SFMO(4,JJK,ISA)+1.
      JRT=1
      RETURN
      END