      SUBROUTINE NITVOL(Z5)
!     APEX0806
!     THIS SUBPROGRAM SIMULATES THE TRANSFORMATION FROM NH3 TO SOL N, 
!     AND THE VOLATILIZATION OF NH3 USING MODIFIED METHODS OF REDDY AND 
!     OF THE CERES MODEL.
      USE PARM
      X1=.41*(STMP(ISL,ISA)-5.)
      IF(X1<=0.)RETURN
      IF(ISL/=LID(1,ISA))THEN
          FCEC=MAX(PRMT(27),1.-.038*CEC(ISL,ISA))
          FZ=1.-Z5/(Z5+EXP(SCRP(12,1)-SCRP(12,2)*Z5))
          AKAV=X1*FCEC*FZ
      ELSE
          FAF=.335+.16*LOG(U10(IRF(ISA))+.2)
          AKAV=X1*FAF
      END IF
      IF(PH(ISL,ISA)>7.)THEN
          IF(PH(ISL,ISA)>7.4)THEN
              FPH=5.367-.599*PH(ISL,ISA)
          ELSE
              FPH=1.
          END IF
      ELSE
          FPH=.307*PH(ISL,ISA)-1.269
      END IF    
      AKAN=X1*SUT*FPH
      AKAV=AKAV*SUT
      XX=AKAV+AKAN
      IF(XX>0.)THEN
          F=MIN(PRMT(80),1.-EXP(-XX))
          X1=F*WNMA(ISL,ISA)
          !X2=1.-Z5/(Z5+EXP(5.-.05*Z5))
          AVOL=X1*PRMT(72)
          RNIT=X1-AVOL
          WNMA(ISL,ISA)=WNMA(ISL,ISA)-AVOL-RNIT
          WNMN(ISL,ISA)=WNMN(ISL,ISA)+RNIT
          SVOL=SVOL+AVOL
          SNIT=SNIT+RNIT
      END IF
      RETURN
      END