      SUBROUTINE NFIX
!     APEX0806
!     THIS SUBPROGRAM ESTIMATES N FIXATION FOR LEGUMES.
      USE PARM
      IF(HUI(JJK,ISA)<.15.OR.HUI(JJK,ISA)>.75)GO TO 8
      SUM=0.
      TOT=0.
      ADD=0.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          IF(Z(ISL,ISA)>.3)GO TO 3
          SUM=SUM+ST(ISL,ISA)-S15(ISL,ISA)
          TOT=TOT+FC(ISL,ISA)-S15(ISL,ISA)
      END DO
      GO TO 4
    3 L1=LID(J-1,ISA)
      RTO=(.3-Z(L1,ISA))/(Z(ISL,ISA)-Z(L1,ISA))
      SUM=SUM+(ST(ISL,ISA)-S15(ISL,ISA))*RTO
      TOT=TOT+(FC(ISL,ISA)-S15(ISL,ISA))*RTO
    4 X1=SUM/TOT
      IF(X1<=.25)GO TO 8
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          IF(Z(ISL,ISA)>RD(JJK,ISA))GO TO 6
          ADD=ADD+WNMN(ISL,ISA)
      END DO
      GO TO 7
    6 L1=LID(J-1,ISA)
      RTO=(RD(JJK,ISA)-Z(L1,ISA))/(Z(ISL,ISA)-Z(L1,ISA))
      ADD=ADD+WNMN(ISL,ISA)*RTO
    7 FXN=1.5-.005*ADD/RD(JJK,ISA)
      IF(FXN<=0.)GO TO 8
      FXW=1.333*X1-.333
      FXG=(HUI(JJK,ISA)-.1)*5.
      FXS=4.-5.*HUI(JJK,ISA)
      FXP=MIN(FXG,FXS,1.)
      FIXR=MIN(FXW,FXN,1.)*FXP
      WFX=FIXR*UNM
    8 WFX=PRMT(7)*WFX+(1.-PRMT(7))*UNM
      WFX=MIN(PRMT(28),WFX)
      SMM(43,MO,ISA)=SMM(43,MO,ISA)+WFX
      UNM=UNM-WFX
      RETURN
      END