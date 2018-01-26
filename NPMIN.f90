      SUBROUTINE NPMIN
!     APEX0806
!     THIS SUBPROGRAM COMPUTES P FLUX BETWEEN THE SOL, ACTIVE MINERAL
!     AND STABLE MINERAL P POOLS.
      USE PARM
      RTO=MIN(.8,PSP(ISL,ISA)/(1.-PSP(ISL,ISA)))
      RMN=PRMT(84)*(WPML(ISL,ISA)-WPMA(ISL,ISA)*RTO)
      X1=4.*WPMA(ISL,ISA)-WPMS(ISL,ISA)
      IF(X1>500.)THEN
          ROC=10.**(LOG10(BK(ISL,ISA))+LOG10(X1))
      ELSE
          ROC=BK(ISL,ISA)*X1
      END IF
      ROC=PRMT(85)*ROC
      WPMS(ISL,ISA)=WPMS(ISL,ISA)+ROC
      WPMA(ISL,ISA)=WPMA(ISL,ISA)-ROC+RMN
      WPML(ISL,ISA)=WPML(ISL,ISA)-RMN
      RETURN
      END