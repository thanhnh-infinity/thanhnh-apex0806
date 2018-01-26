      SUBROUTINE HSWU(CPWU,RGS)
!     APEX0806
!     THIS SUBPROGRAM DISTRIBUTES PLANT EVAPORATION THROUGH THE ROOT
!     ZONE AND CALCULATES ACTUAL PLANT WATER USE BASED ON SOIL WATER
!     AVAILABILITY.
      USE PARM
      BLM=S15(ISL,ISA)
      IF(Z(ISL,ISA)<=.5)BLM=PRMT(5)*S15(ISL,ISA)
      IF(ISL/=LID(1,ISA))THEN
          CALL CRGBD(RGS)
          CPWU=CPWU*RGS
      END IF
      SUM=EP(JJK)*(1.-EXP(-UB1(ISA)*GX/RD(JJK,ISA)))/UOB(ISA)
      TOS=36.*ECND(ISL,ISA)
      XX=LOG10(S15(ISL,ISA))
      X1=3.1761-1.6576*(LOG10(ST(ISL,ISA))-XX)/(LOG10(FC(ISL,ISA))-XX)
      IF(X1<4.)THEN
          WTN=MAX(5.,10.**X1)
      ELSE
          WTN=10000.
          GO TO 4
      END IF
      XX=TOS+WTN
      IF(XX>5000.)GO TO 4
      F=1.-XX/(XX+EXP(SCRP(22,1)-SCRP(22,2)*XX))
      UW(ISL)=MIN(SUM-CPWU*AEP(JJK)-(1.-CPWU)*UX,ST(ISL,ISA)-BLM)*F*RGS
      UW(ISL)=MAX(0.,UW(ISL))
    4 UX=SUM
      RETURN
      END