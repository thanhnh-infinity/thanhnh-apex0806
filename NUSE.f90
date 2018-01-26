      SUBROUTINE NUSE
!     APEX0806
!     THIS SUBPROGRAM CALCULATES THE DAILY POTENTIAL SOIL SUPPLY OF P
!     FOR EACH LAYER.
      USE PARM
      XX=1.5*UPM/RW(JJK,ISA)
      DO J=1,LRD(ISA)
          ISL=LID(J,ISA)
          UN(ISL)=WNMN(ISL,ISA)*UW(ISL)/(ST(ISL,ISA)+.001)
          SUN=SUN+UN(ISL)
          F=1000.*WPML(ISL,ISA)/WT(ISL,ISA)
          IF(F>30.)THEN
              F=1.
          ELSE
              F=F/(F+EXP(SCRP(11,1)-SCRP(11,2)*F))
          END IF
          UP(ISL)=XX*F*RWT(ISL,JJK,ISA)
          IF(UP(ISL)>=WPML(ISL,ISA))UP(ISL)=.9*WPML(ISL,ISA)
          SUP=SUP+UP(ISL)
      END DO
      RETURN
      END