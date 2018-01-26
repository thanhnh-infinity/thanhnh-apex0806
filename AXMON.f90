      SUBROUTINE AXMON
!     APEX0806
!     THIS SUBPROGRAM DETERMINES THE MONTH, GIVEN THE DAY OF THE YEAR.
      USE PARM
      IF(JDA>NC(2))THEN
          M=MO
          DO MO=M,12
              M1=MO+1
              NDA=NC(M1)-NYD
              IF(JDA<=NDA)RETURN
          END DO
      END IF
          MO=1
      RETURN
      END