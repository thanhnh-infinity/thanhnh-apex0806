      SUBROUTINE HRFDTS
      ! THIS SUBPROGRAM GENERATES RAINFALL DISTRIBUTIONS AT DTHY TIME 
      ! INTERVALS USING AN S CURVE
      USE PARM
      T1=0.
      DTX=DTHY/.12
      X1=RFV(IRF(ISA))
      I=1
      DO 
          T1=T1+DTX
          RFF=T1/(T1+EXP(SCRP(26,1)-SCRP(26,2)*T1))
          I=I+1
          RFDT(I)=X1*RFF
          IF(T1>100.)EXIT
      END DO
      NRF=I     
      RTO=X1/RFDT(NRF)
      DO J=1,NRF
          RFDT(J)=RFDT(J)*RTO
      END DO
      RETURN
      END