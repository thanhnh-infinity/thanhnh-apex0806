      SUBROUTINE NUP
!     APEX0806
!     THIS SUBPROGRAM CALCULATES THE DAILY N DEMAND FOR OPTIMAL PLANT
!     GROWTH.
      USE PARM
      IF(NUPC/=0)THEN
          CNT=(BN(4,JJK)-BN(3,JJK))*(1.-HUI(JJK,ISA)/(HUI(JJK,ISA)+EXP(&
          &BN(1,JJK)-BN(2,JJK)*HUI(JJK,ISA))))+BN(3,JJK)
      ELSE
          CNT=BN(2,JJK)+BN(1,JJK)*EXP(-BN(4,JJK)*HUI(JJK,ISA))
      END IF
      UN2=CNT*DM(JJK,ISA)*1000.
      IF(UN2<UN1(JJK,ISA))UN2=UN1(JJK,ISA)
      UNM=MIN(4000.*BN(3,JJK)*DDM(JJK),UN2-UN1(JJK,ISA))
      UNM=MAX(0.,UNM)
      RETURN
      END