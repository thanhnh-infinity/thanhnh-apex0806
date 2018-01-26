      SUBROUTINE NUTS(U1,U2,UU)
!     APEX0806
!     THIS SUBPROGRAM CALCULATES THE PLANT STRESS FACTOR CAUSED BY LIMIT
!     SUPPLY OF N OR P.
      USE PARM
      UU=200.*(U1/(U2+1.E-10)-.5)
      IF(UU>0.)THEN
	      IF(UU>200.)THEN
	          UU=1.
          ELSE
              UU=UU/(UU+EXP(SCRP(8,1)-SCRP(8,2)*UU))
	      END IF
      ELSE
          UU=0.
      END IF    
      RETURN
      END