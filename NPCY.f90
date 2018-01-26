      SUBROUTINE NPCY(SN3I)
!     APEX0806
!     THIS SUBPROGRAM IS THE MASTER NUTRIENT CYCLING SUBROUTINE.
!     CALLS NPMIN, NYNIT, NLCH, NCNMI, AND NDNIT FOR EACH SOIL
!     LAYER.
      USE PARM 
      STDX=0.
      DO K=1,LC
          STDX=STDX+STD(K,ISA)
      END DO          
      SGMN=0.
      SDN=0.
      SN2=0.
      SMP=0.
!     CALL NCONT
      TSFN(IDO)=0.
      QRFN(IDO)=0.
      SVOL=0.
      SNIT=0.
      TRSP=0.
      TSFS=0.
      XX=0.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          RSPC(ISL)=0.
          RNMN(ISL)=0.
          WDN=0.
          X1=ST(ISL,ISA)-S15(ISL,ISA)
          IF(X1<0.)THEN
              SUT=.1*(ST(ISL,ISA)/S15(ISL,ISA))**2
          ELSE
              SUT=.1+.9*SQRT(X1/(FC(ISL,ISA)-S15(ISL,ISA)))
          END IF    
          CALL NPMIN
          IF(ISL/=LID(1,ISA))THEN
              CALL NLCH(SN3I)
          ELSE
              ZZ=MIN(.9,PRMT(76)*(1.+.1*RFV(IRF(ISA))))
              IF(STDX+STDO(ISA)>.001)CALL NSTDFAL(ZZ)
              CALL NYNIT(SN3I)
          END IF    
          IF(ISL/=IDR(ISA))THEN
              QRFN(IDO)=QRFN(IDO)+QSFN
              SMM(84,MO,ISA)=SMM(84,MO,ISA)+QSFN
          ELSE
              SMM(47,MO,ISA)=SMM(47,MO,ISA)+QSFN
              QDRN(IDO)=QSFN
          END IF
          TSFN(IDO)=TSFN(IDO)+SSFN
          TSFS=TSFS+SSST
          Z5=500.*(Z(ISL,ISA)+XX)
          IF(WNMA(ISL,ISA)>.01)CALL NITVOL(Z5)
          IF(STMP(ISL,ISA)>0.)THEN
              CDG=STMP(ISL,ISA)/(STMP(ISL,ISA)+EXP(SCRP(14,1)-SCRP(14,2)*&
              STMP(ISL,ISA)))
              IF(RZ(ISA)>=XX)THEN
                  CALL NCNMI(Z5,CS)
                  CALL NPMN(CS)
                  SMP=SMP+WMP
              END IF
              DZ=1000.*(Z(ISL,ISA)-XX)
              IF(IDNT>0)THEN
                  CALL NDNITAK(DZ)
              ELSE
                  IF(ST(ISL,ISA)>FC(ISL,ISA))CALL NDNIT
              END IF
              SDN=SDN+WDN
              SN2=SN2+DN2O                  
          END IF
          XX=Z(ISL,ISA)
      END DO
      GWSN(ISA)=GWSN(ISA)+SSO3(LNS)
      SMM(80,MO,ISA)=SMM(80,MO,ISA)+RSFN(IDO)
      SMM(51,MO,ISA)=SMM(51,MO,ISA)+VAP(ISA)+VPU(ISA)
      VAR(51,ISA)=VAP(ISA)+VPU(ISA)
      SMM(133,MO,ISA)=SMM(133,MO,ISA)+VSLT(ISA) 
      RETURN
      END