      FUNCTION CAHU(J,K,BASE,NHS)
!     APEX0806
!     THIS SUBPROGRAM ACCUMULATES HEAT UNITS FOR USE IN CPTHU.
      USE PARM
      CAHU=0.
      MO=1
      DO JDA=J,K
          CALL AXMON
          IF(JDHU<=366)THEN
              CALL WHLRMX(JDA)
              IF(HRLT<WDRM(ISA).AND.NHS==0)CYCLE
          END IF
          TA=ARALT(TAV,XX)
          TGX=TA-BASE
          IF(TGX>0.)CAHU=CAHU+TGX
      END DO
      RETURN
      END