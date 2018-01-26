      SUBROUTINE WSOLRA(I)
!     APEX0806
!     THIS SUBPROGRAM SIMULATES DAILY SOLAR RADIATION FROM A NORMAL
!     DISTRIBUTION.
      USE PARM 
      RX=RAMX-RM
      SRAD(I)=RM+WX(3)*RX/4.
      IF(SRAD(I)<=0.)SRAD(I)=.05*RAMX
      RETURN
      END
