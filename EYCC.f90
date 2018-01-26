      SUBROUTINE EYCC
!     APEX0806
!     THIS SUBPROGRAM ESTIMATES THE USLE C FACTOR BASED ON PLANT POP &
!     BIOMASS & RESIDUE COVER
      USE PARM
      NN=NCP(IRO(ISA),ISA)
      SUM=0.
      TOT=0.
      DO I=1,NBSL(ISA)
          L=LID(I,ISA)
          IF(Z(L,ISA)>PMX(ISA))GO TO 3
          IF(I>1)TOT=TOT+RSD(L,ISA)
          DO K=1,NN
              SUM=SUM+RWT(L,JE(K,ISA),ISA)
          END DO
      END DO
      GO TO 4
    3 KK=LID(I-1,ISA)
      RTO=(PMX(ISA)-Z(KK,ISA))/(Z(L,ISA)-Z(KK,ISA))
      TOT=TOT+RTO*RSD(L,ISA)
      DO K=1,NN
          SUM=SUM+RTO*RWT(L,JE(K,ISA),ISA)
      END DO
    4 SUM=SUM/PMX(ISA)
      TOT=TOT/(PMX(ISA)-.01)
      FRUF=MIN(1.,EXP(-.026*(RRUF(ISA)-6.1)))
      IF(CVRS(ISA)<15.)THEN
          FRSD=EXP(-PRMT(46)*CVRS(ISA))
      ELSE
          FRSD=.0001
      END IF    
      X1=MAX(FGC(ISA),FGSL(ISA))
      FBIO=1.-X1*EXP(-PRMT(47)*CPHT(JJK,ISA))
      FPPL=.9*(1.-CVP(ISA))+.1
      CVF(ISA)=FRSD*FBIO*FRUF*FPPL
      CVF(ISA)=CVF(ISA)*EXP(-.03*ROK(LID(1,ISA),ISA))
      RETURN
      END