      SUBROUTINE CAGRO
!     APEX0806
!     THIS SUBPROGRAM CALCULATES THE DAILY INCREASE IN PLANT BIOMASS,
!     ROOT WEIGHT, AND YIELD BY ADJUSTING THE POTENTIAL VALUES WITH THE
!     ACTIVE STRESS CONSTRAINT.
      USE PARM
      LD1=LID(1,ISA)
      X1=REG(JJK,ISA)*AJWA*SHRL
      RWL=RW(JJK,ISA)
      RGD=DDM(JJK)*X1
      X2=100.*HUI(JJK,ISA)
      AJHI(JJK,ISA)=HI(JJK)*X2/(X2+EXP(SCRP(3,1)-SCRP(3,2)*X2))
      X3=DM(JJK,ISA)-DDM(JJK)
      X4=MAX(1.E-5,X3+RGD)
      DM(JJK,ISA)=X4
      DM1(JJK,ISA)=DM1(JJK,ISA)+RGD
      X1=HUI(JJK,ISA)
      RF=MAX(.2,RWPC(1,JJK)*(1.-X1)+RWPC(2,JJK)*X1)
      RF=MIN(RF,.99)
	  RW(JJK,ISA)=RF*DM(JJK,ISA)
	  DRW=RW(JJK,ISA)-RWL
      STL(JJK,ISA)=DM(JJK,ISA)-RW(JJK,ISA)
      IF(IDC(JJK)==NDC(7).OR.IDC(JJK)==NDC(8).OR.IDC(JJK)==NDC(10))THEN
          X1=MO*MO
	      FALF=STL(JJK,ISA)*X1*XMTU(JJK)
          SMM(109,MO,ISA)=SMM(109,MO,ISA)+FALF
!         SLAI(JJK,ISA)=MAX(.05,SLAI(JJK,ISA)-FALF)
          DM(JJK,ISA)=DM(JJK,ISA)-FALF
          STL(JJK,ISA)=STL(JJK,ISA)-FALF
          X11=FALF*1000.
          X5=CNY(JJK)*X11
          X6=CPY(JJK)*X11
          IF(FALF>1.E-5)CALL NCNSTD(.1,FALF,X5,LD1)
          FOP(LD1,ISA)=FOP(LD1,ISA)+X6
          UN1(JJK,ISA)=UN1(JJK,ISA)-X5
          UP1(JJK,ISA)=UP1(JJK,ISA)-X6
      END IF
      IF(IDC(JJK)==NDC(3).OR.IDC(JJK)==NDC(6))THEN
          X7=.01*(HUI(JJK,ISA)+.01)**10*STL(JJK,ISA)
          STL(JJK,ISA)=STL(JJK,ISA)-X7
          DM(JJK,ISA)=DM(JJK,ISA)-X7
          STD(JJK,ISA)=STD(JJK,ISA)+X7
          STDL(JJK,ISA)=STDL(JJK,ISA)+CLG(ISA)*X7
          X8=X7*BN(3,JJK)
          XUN=UN1(JJK,ISA)
          IF(XUN-X8<.01)X8=XUN-.01
          X9=X7*BP(3,JJK)
          XUP=UP1(JJK,ISA)
          IF(XUP-X9<.01)X9=XUP-.01
          STDN(JJK,ISA)=STDN(JJK,ISA)+X8
          STDP(JJK,ISA)=STDP(JJK,ISA)+X9
          UN1(JJK,ISA)=XUN-X8
          UP1(JJK,ISA)=XUP-X9
          IF(HUI(JJK,ISA)>.6.AND.STL(JJK,ISA)<.1)HU(JJK,ISA)=0.
      END IF
      SUM=0.
      DO J=1,LRD(ISA)
          ISL=LID(J,ISA)
          RTO=WNMU(ISL,ISA)/(WNMU(ISL,ISA)+WNMN(ISL,ISA)+1.E-5)
          UU=UN(ISL)*RTO
          WNMN(ISL,ISA)=MAX(0.,WNMN(ISL,ISA)-UN(ISL)+UU)
          WNMU(ISL,ISA)=MAX(0.,WNMU(ISL,ISA)-UU)
          IF(DRW<0.)THEN      
              UTO=RWT(ISL,JJK,ISA)/RWL
              X1=-DRW*UTO
              X2=MIN(UN1(JJK,ISA),1000.*BN(3,JJK)*X1)
              CALL NCNSTD(.1,X1,X2,ISL)
              UN1(JJK,ISA)=UN1(JJK,ISA)-X2
          ELSE
              UTO=UW(ISL)/(AEP(JJK)+1.E-20)
          END IF
          ST(ISL,ISA)=MAX(1.E-10,ST(ISL,ISA)-UW(ISL))
          IF(WPML(ISL,ISA)>UP(ISL))THEN
              WPML(ISL,ISA)=WPML(ISL,ISA)-UP(ISL)
          ELSE
              UP(ISL)=WPML(ISL,ISA)
              WPML(ISL,ISA)=0.
          END IF
          RWT(ISL,JJK,ISA)=RWT(ISL,JJK,ISA)+DRW*UTO
          SUM=SUM+RWT(ISL,JJK,ISA)
      END DO
      RW(JJK,ISA)=SUM
      RETURN
      END