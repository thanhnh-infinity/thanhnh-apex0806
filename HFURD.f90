      SUBROUTINE HFURD
!     APEX0806
!     THIS SUBPROGRAM COMPUTES THE STORAGE VOLUME OF FURROW DIKES GIVEN
!     DIKE INTERVAL AND HEIGHT, RIDGE HEIGHT, AND SLOPE.
      USE PARM
      S=STP(ISA)
      DH=DHT(ISA)*.001
      H2=2.*DH
      X1=.001*DKHL(ISA)
      TW=RINT(ISA)-X1
      BW=MAX(TW-4.*X1,.1*TW)
      DQI=DKIN(ISA)-X1
      D2=DH*(1.-2.*S)
      D3=DH-S*(DQI-H2)
      X1=(TW-BW)/DH
      TW2=BW+D2*X1
      TW3=BW+D3*X1
      A2=.5*D2*(TW2+BW)
      A3=.5*D3*(TW3+BW)
      XX=DH/S
      ZZ=DQI-H2
      X1=FDSF(ISA)/(RINT(ISA)*DKIN(ISA))
      IF(XX>ZZ)THEN
          DVOL=X1*(A2*DH+.5*(A2+A3)*(DQI-4.*DH)+A3*D3)
      ELSE
          DVOL=X1*A2*(DH+.5*(XX-H2))
      END IF
      RETURN
      END