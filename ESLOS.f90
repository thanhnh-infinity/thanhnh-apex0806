      SUBROUTINE ESLOS(JRT)
!     APEX0806
!     THIS SUBPROGRAM CALCULATES THE THICKNESS OF SOIL REMOVED BY EACH
!     EROSION EVENT, MOVES THE TOP LAYER INTO THE SECOND LAYER BY A DIST
!     EQUAL THE ERODED THICKNESS, AND ADJUSTS THE TOP LAYER SOIL
!     PROPERTIES BY INTERPOLATION.  WHEN THE SECOND LAYER IS REDUCED TO
!     TO A THICKNESS OF 10 MM, IT IS PLACED INTO THE THIRD LAYER.
      USE PARM
	  WFN(X1,X2,Z1,Z2,ZZ)=(Z1*X1+Z2*X2)/ZZ
	  JRT=0
      TK=.1*YERO/BD(LID(1,ISA),ISA)
      THK(ISA)=THK(ISA)+TK
      TK=TK*.001
      X1=TK+TKR(ISA)
      IF(X1>.01)THEN
          TKR(ISA)=X1-.01
          TK=.01
      ELSE
          TKR(ISA)=0.
      END IF
      W1=Z(LID(1,ISA),ISA)-TK
      W2=TK
      WT(LID(1,ISA),ISA)=WT(LID(1,ISA),ISA)-YERO
      IF(Z(LID(2,ISA),ISA)-Z(LID(1,ISA),ISA)-TK>.01)GO TO 10
!     REMOVE LAYER 2 AND PLACE SMALL REMAINING CONTENTS IN LAYER 3
      ISL=LID(3,ISA)
      K=LID(2,ISA)
      Z1=Z(K,ISA)-Z(LID(1,ISA),ISA)
      Z2=Z(ISL,ISA)-Z(K,ISA)
      VV=Z1+Z2
	  ZZ=.01/VV
      PSP(ISL,ISA)=WFN(PSP(K,ISA),PSP(ISL,ISA),Z1,Z2,VV)
!     BDM(ISL,ISA)=WFN(BDM(K,ISA),BDM(ISL,ISA),Z1,Z2,VV)
      CLA(ISL,ISA)=WFN(CLA(K,ISA),CLA(ISL,ISA),Z1,Z2,VV)
      SIL(ISL,ISA)=WFN(SIL(K,ISA),SIL(ISL,ISA),Z1,Z2,VV)
      SAN(ISL,ISA)=100.-CLA(ISL,ISA)-SIL(ISL,ISA)
      ROK(ISL,ISA)=WFN(ROK(K,ISA),ROK(ISL,ISA),Z1,Z2,VV)
      SATC(ISL,ISA)=WFN(SATC(K,ISA),SATC(ISL,ISA),Z1,Z2,VV)
      PH(ISL,ISA)=WFN(PH(K,ISA),PH(ISL,ISA),Z1,Z2,VV)
      WT(ISL,ISA)=WT(ISL,ISA)+WT(K,ISA)
      BD(ISL,ISA)=.0001*WT(ISL,ISA)/VV
      WCMU(ISL,ISA)=WCMU(ISL,ISA)+WCMU(K,ISA)
      WCOU(ISL,ISA)=WCOU(ISL,ISA)+WCOU(K,ISA)
      WNMN(ISL,ISA)=WNMN(ISL,ISA)+WNMN(K,ISA)
      WNMA(ISL,ISA)=WNMA(ISL,ISA)+WNMA(K,ISA)
      WNMU(ISL,ISA)=WNMU(ISL,ISA)+WNMU(K,ISA)
	  WHPN(ISL,ISA)=WHPN(ISL,ISA)+WHPN(K,ISA)
      WHSN(ISL,ISA)=WHSN(ISL,ISA)+WHSN(K,ISA)
      WBMN(ISL,ISA)=WBMN(ISL,ISA)+WBMN(K,ISA)
      WLSN(ISL,ISA)=WLSN(ISL,ISA)+WLSN(K,ISA)
      WLMN(ISL,ISA)=WLMN(ISL,ISA)+WLMN(K,ISA)
      WHPC(ISL,ISA)=WHPC(ISL,ISA)+WHPC(K,ISA)
      WHSC(ISL,ISA)=WHSC(ISL,ISA)+WHSC(K,ISA)
      WBMC(ISL,ISA)=WBMC(ISL,ISA)+WBMC(K,ISA)
      WLS(ISL,ISA)=WLS(ISL,ISA)+WLS(K,ISA)
      WLM(ISL,ISA)=WLM(ISL,ISA)+WLM(K,ISA)
      WLSL(ISL,ISA)=WLSL(ISL,ISA)+WLSL(K,ISA)
      WLSC(ISL,ISA)=WLSC(ISL,ISA)+WLSC(K,ISA)
      WLMC(ISL,ISA)=WLMC(ISL,ISA)+WLMC(K,ISA)
      WLSLC(ISL,ISA)=WLSLC(ISL,ISA)+WLSLC(K,ISA)
      WLSLNC(ISL,ISA)=WLSC(ISL,ISA)-WLSLC(ISL,ISA)
      WPO(ISL,ISA)=WPO(ISL,ISA)+WPO(K,ISA)
      RSDM(ISL,ISA)=RSDM(ISL,ISA)+RSDM(K,ISA)
      WNOU(ISL,ISA)=WNOU(ISL,ISA)+WNOU(K,ISA)
      WPML(ISL,ISA)=WPML(ISL,ISA)+WPML(K,ISA)
	  WPMA(ISL,ISA)=WPMA(ISL,ISA)+WPMA(K,ISA)
      WPMU(ISL,ISA)=WPMU(ISL,ISA)+WPMU(K,ISA)
	  FOP(ISL,ISA)=FOP(ISL,ISA)+FOP(K,ISA)
      WPOU(ISL,ISA)=WPOU(ISL,ISA)+WPOU(K,ISA)
      WPMS(ISL,ISA)=WPMS(ISL,ISA)+WPMS(K,ISA)
      S15(ISL,ISA)=S15(ISL,ISA)+S15(K,ISA)
      FC(ISL,ISA)=FC(ISL,ISA)+FC(K,ISA)
      PO(ISL,ISA)=PO(ISL,ISA)+PO(K,ISA)
      CALL SPOFC(ISL)
      ST(ISL,ISA)=ST(ISL,ISA)+ST(K,ISA)
      LORG(LID(1,ISA),ISA)=LORG(K,ISA)
!     SPLIT LAYER NEAREST SURFACE WITH THICKNESS > 0.15 M IN HALF
      ZMX=0.
      L1=LID(1,ISA)
      MXZ=K
      DO J=3,NBSL(ISA)
          L=LID(J,ISA)
          ZZ=Z(L,ISA)-Z(L1,ISA)
          IF(ZZ>=.15)THEN
              MXZ=J
              ZMX=ZZ
              GO TO 7
          ELSE
              L1=L
              IF(ZZ<=ZMX+.01)CYCLE
              ZMX=ZZ
              MXZ=J
          END IF
      END DO
      L=LID(MXZ,ISA)
      L1=LID(MXZ-1,ISA)
      IF(ZMX>ZQT)GO TO 7
      NBSL(ISA)=NBSL(ISA)-1
      IF(NBSL(ISA)>2)GO TO 1
	  JRT=1
      RETURN 
    1 XNS(ISA)=NBSL(ISA)-1
      DO ISL=2,NBSL(ISA)
          LID(ISL,ISA)=LID(ISL+1,ISA)
      END DO
      GO TO 10
    7 MX1=MXZ-1
      IF(MX1/=2)THEN
          DO ISL=2,MX1
              LID(ISL,ISA)=LID(ISL+1,ISA)
          END DO
      END IF
      LID(MX1,ISA)=K
      LORG(K,ISA)=LORG(LID(MXZ,ISA),ISA)
      CALL SPLA(L,L1,K,.5,0)
!     ADJUST LAYER 1 BY INTERPOLATING BETWEEN LAYER 2 USING ERODED THICK
   10 ISL=LID(2,ISA)
      L1=LID(1,ISA)
      PSP(L1,ISA)=WFN(PSP(L1,ISA),PSP(ISL,ISA),W1,W2,Z(L1,ISA))
!     BDM(L1,ISA)=WFN(BDM(L1,ISA),BDM(ISL,ISA),W1,W2,Z(L1,ISA))
      CLA(L1,ISA)=WFN(CLA(L1,ISA),CLA(ISL,ISA),W1,W2,Z(L1,ISA))
      SIL(L1,ISA)=WFN(SIL(L1,ISA),SIL(ISL,ISA),W1,W2,Z(L1,ISA))
      SAN(L1,ISA)=100.-CLA(L1,ISA)-SIL(L1,ISA)
      BD(L1,ISA)=WFN(BD(L1,ISA),BD(ISL,ISA),W1,W2,Z(L1,ISA))
      SATC(L1,ISA)=WFN(SATC(L1,ISA),SATC(ISL,ISA),W1,W2,Z(L1,ISA))
      PH(L1,ISA)=WFN(PH(L1,ISA),PH(ISL,ISA),W1,W2,Z(L1,ISA))
      XX=Z(ISL,ISA)-Z(L1,ISA)
      GX=XX-TK
      RTO=GX/XX
      RX=1.
      IF(ROK(L1,ISA)>1.E-10)THEN
          X3=ROK(L1,ISA)
          !ROK(L1,ISA)=.01*X3/W1
          ROK(L1,ISA)=WFN(ROK(L1,ISA),ROK(ISL,ISA),W1,W2,Z(L1,ISA))
          RX=(100.-ROK(L1,ISA))/(100.-X3)
      END IF
      W1=W1/Z(L1,ISA)
      W2=W2/XX
      S15(L1,ISA)=WFN(S15(L1,ISA),S15(ISL,ISA),W1,W2,1.)*RX
      FC(L1,ISA)=WFN(FC(L1,ISA),FC(ISL,ISA),W1,W2,1.)*RX
      PO(ISL,ISA)=PO(ISL,ISA)*RTO
      S15(ISL,ISA)=S15(ISL,ISA)*RTO
      FC(ISL,ISA)=FC(ISL,ISA)*RTO
      CALL SPOFC(ISL)
      WT(ISL,ISA)=WT(ISL,ISA)*RTO
      PO(L1,ISA)=WFN(PO(L1,ISA),PO(ISL,ISA),W1,W2,1.)*RX
      CALL SPOFC(L1)
      WCMU(L1,ISA)=WCMU(L1,ISA)+EAJL(WCMU(ISL,ISA),W2)
      WCOU(L1,ISA)=WCOU(L1,ISA)+EAJL(WCOU(ISL,ISA),W2)
      WNMN(L1,ISA)=WNMN(L1,ISA)+EAJL(WNMN(ISL,ISA),W2)
      WNMA(L1,ISA)=WNMA(L1,ISA)+EAJL(WNMA(ISL,ISA),W2)
      WNMU(L1,ISA)=WNMU(L1,ISA)+EAJL(WNMU(ISL,ISA),W2)
      WPO(L1,ISA)=WPO(L1,ISA)+EAJL(WPO(ISL,ISA),W2)
      RSDM(L1,ISA)=RSDM(L1,ISA)+EAJL(RSDM(ISL,ISA),W2)
      WNOU(L1,ISA)=WNOU(L1,ISA)+EAJL(WNOU(ISL,ISA),W2)
      WPML(L1,ISA)=WPML(L1,ISA)+EAJL(WPML(ISL,ISA),W2)
	  WPMA(L1,ISA)=WPMA(L1,ISA)+EAJL(WPMA(ISL,ISA),W2)
      WPMU(L1,ISA)=WPMU(L1,ISA)+EAJL(WPMU(ISL,ISA),W2)
	  WHPN(L1,ISA)=WHPN(L1,ISA)+EAJL(WHPN(ISL,ISA),W2)
      WHSN(L1,ISA)=WHSN(L1,ISA)+EAJL(WHSN(ISL,ISA),W2)
      WBMN(L1,ISA)=WBMN(L1,ISA)+EAJL(WBMN(ISL,ISA),W2)
      WLSN(L1,ISA)=WLSN(L1,ISA)+EAJL(WLSN(ISL,ISA),W2)
      WLMN(L1,ISA)=WLMN(L1,ISA)+EAJL(WLMN(ISL,ISA),W2)
      WHPC(L1,ISA)=WHPC(L1,ISA)+EAJL(WHPC(ISL,ISA),W2)
      WHSC(L1,ISA)=WHSC(L1,ISA)+EAJL(WHSC(ISL,ISA),W2)
      WBMC(L1,ISA)=WBMC(L1,ISA)+EAJL(WBMC(ISL,ISA),W2)
      WLS(L1,ISA)=WLS(L1,ISA)+EAJL(WLS(ISL,ISA),W2)
      WLM(L1,ISA)=WLM(L1,ISA)+EAJL(WLM(ISL,ISA),W2)
      WLSL(L1,ISA)=WLSL(L1,ISA)+EAJL(WLSL(ISL,ISA),W2)
      WLSC(L1,ISA)=WLSC(L1,ISA)+EAJL(WLSC(ISL,ISA),W2)
      WLMC(L1,ISA)=WLMC(L1,ISA)+EAJL(WLMC(ISL,ISA),W2)
      WLSLC(L1,ISA)=WLSLC(L1,ISA)+EAJL(WLSLC(ISL,ISA),W2)
      WLSLNC(L1,ISA)=WLSC(L1,ISA)-WLSLC(L1,ISA)
	  FOP(L1,ISA)=FOP(L1,ISA)+EAJL(FOP(ISL,ISA),W2)
      WPOU(L1,ISA)=WPOU(L1,ISA)+EAJL(WPOU(ISL,ISA),W2)
      WPMS(L1,ISA)=WPMS(L1,ISA)+EAJL(WPMS(ISL,ISA),W2)
      ST(L1,ISA)=ST(L1,ISA)+EAJL(ST(ISL,ISA),W2)
      XY=.01*(1.-.01*ROK(L1,ISA))
      XX=XY
      WT(L1,ISA)=BD(L1,ISA)*100.
      ABD(ISA)=WT(L1,ISA)
      DO I=2,NBSL(ISA)
          ISL=LID(I,ISA)
          Z(ISL,ISA)=Z(ISL,ISA)-TK
          ABD(ISA)=ABD(ISA)+WT(ISL,ISA)
      END DO
      XX=Z(LID(NBSL(ISA),ISA),ISA)
      IF(XX<=ZF)THEN
	      JRT=1
          RETURN 
      ELSE
          IF(BIG(ISA)>=XX)THEN
              BIG(ISA)=XX
              PMX(ISA)=XX
          END IF
          ABD(ISA)=ABD(ISA)*1.E-4/Z(LID(NBSL(ISA),ISA),ISA)
      END IF
      RETURN
      END