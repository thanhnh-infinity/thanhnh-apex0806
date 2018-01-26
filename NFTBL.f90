      SUBROUTINE NFTBL(K)
!     APEX0806
!     THIS SUBPROGRAM READS FERTILIZER TABLE TO DETERMINE PARAMETERS OF
!     INPUT FERTILIZER
      USE PARM
      DIMENSION YTP(10)
      IF(NDF>0)THEN
          DO K=1,NDF
              IF(KDF(K)==JX(7))RETURN
          END DO
      END IF
      NDF=NDF+1
      KDF(NDF)=JX(7)
      KDF1(JX(7))=NDF
!     READ FERTILIZER TABLE
!  1  FTNM = FERTILIZER NAME
!  2  FN   = MINERAL N FRACTION
!  3  FP   = MINERAL P FRACTION
!  4  FK   = MINERAL K FRACTION
!  5  FNO  = ORGANIC N FRACTION
!  6  FPO  = ORGANIC P FRACTION
!  7  FNH3 = AMMONIA N FRACTION(FNH3/FN)
!  8  FOC  = ORGANIC C FRACTION
!  9  FSLT = SALT FRACTION
! 10  FCST = COST OF FERTILIZER($/KG)
      I=-1
      DO WHILE(I/=JX(7))
          READ(KR(9),395,IOSTAT=NFL)I,FTNM(NDF),(YTP(J),J=2,10)
          IF(NFL/=0)THEN
              WRITE(*,*)'FERT NO = ',JX(7),' NOT IN FERT LIST FILE     &
              SAID = ',NBSA(ISA)
              PAUSE
              STOP
          END IF
      END DO
      FN(NDF)=YTP(2)
      FP(NDF)=YTP(3)
!     FK(NDF)=YTP(4)
      FNO(NDF)=YTP(5)
      FPO(NDF)=YTP(6)
      FNMA(NDF)=YTP(7)
      FOC(NDF)=YTP(8)
      FSLT(NDF)=YTP(9)
      FCST(NDF)=YTP(10)
      REWIND KR(9)
      RETURN
  395 FORMAT(1X,I4,1X,A8,10F8.3)
      END