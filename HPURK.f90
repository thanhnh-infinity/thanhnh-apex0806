      SUBROUTINE HPURK
!     APEX0806
!     THIS SUBPROGRAM IS THE MASTER PERCOLATION COMPONENT.  IT MANAGES
!     THE ROUTING PROCESS
      USE PARM
      ADD=0.
      SUM=0.
      TOT=0.
      XX=0.
      CPVV=0.
      QDR(IDO)=0.
      VAR(98,ISA)=FPF(ISA)
      SEP=RFV(IRF(ISA))-QVOL(IDO)+FPF(ISA)
      LD1=LID(1,ISA)
	  STLT(ISA)=PRMT(51)*RSD(LD1,ISA)
      SWLT(ISA)=SWLT(ISA)+SEP
      IF(SWLT(ISA)>STLT(ISA))THEN
          SEP=SWLT(ISA)-STLT(ISA)
          SWLT(ISA)=STLT(ISA)
      ELSE
          SEP=0.
      END IF
      DO KK=1,NBSL(ISA)
          ISL=LID(KK,ISA)
          DZ=Z(ISL,ISA)-XX
          XX=Z(ISL,ISA)
          ST(ISL,ISA)=ST(ISL,ISA)+SEP
          IF(WTBL(ISA)<=Z(ISL,ISA))THEN
              SSF(ISL,ISA)=0.
              PKRZ(ISL)=0.
              CPFH(ISL,ISA)=0.
              SEP=0.
              CYCLE
          END IF
          CPVV=SEP*CPRV(ISL,ISA)
          X1=SEP-CPVV
          CPVH(IDO)=X1*CPRH(ISL,ISA)
          ST(ISL,ISA)=MAX(1.E-5,ST(ISL,ISA)-CPVV-CPVH(IDO))
	      IF(RSAE(ISA)>0.)THEN
	          SATX=MAX(1.E-10,RSHC(ISA))
	      ELSE
	          SATX=SATC(ISL,ISA)
          END IF
          CALL HPERC(DZ,SATX)
          ST(ISL,ISA)=MAX(1.E-5,ST(ISL,ISA)-SEP-SST(IDO)-QRF(IDO))
          IF(ISL/=IDR(ISA))THEN
              SUM=SUM+QRF(IDO)
          ELSE
              SMM(17,MO,ISA)=SMM(17,MO,ISA)+QRF(IDO)
              VAR(17,ISA)=QRF(IDO)
              QDR(IDO)=QRF(IDO)
          END IF
          ADD=ADD+SST(IDO)
          TOT=TOT+CPVH(IDO)
          SSF(ISL,ISA)=SST(IDO)
          QSF(ISL,ISA)=QRF(IDO)
          CPFH(ISL,ISA)=CPVH(IDO)
          SEP=SEP+CPVV
          PKRZ(ISL)=SEP
      END DO
      SST(IDO)=ADD
      QRF(IDO)=SUM
      CPVH(IDO)=TOT
      L1=LD1
      DO K=NBSL(ISA),2,-1
          ISL=LID(K,ISA)
          L1=LID(K-1,ISA)
          XX=ST(ISL,ISA)-PO(ISL,ISA)
          IF(XX>0.)THEN
              ST(L1,ISA)=ST(L1,ISA)+XX
              PKRZ(L1)=MAX(0.,PKRZ(L1)-XX)
              ST(ISL,ISA)=PO(ISL,ISA)
          END IF
          XX=ST(ISL,ISA)-FC(ISL,ISA)
          IF(XX<=0.)CYCLE
          WP1=LOG10(S15(L1,ISA))
          FC1=LOG10(FC(L1,ISA))
	      IF(ST(L1,ISA)>.01)THEN
	          T1=10.**(3.1761-1.6576*((LOG10(ST(L1,ISA))-WP1)/(FC1-WP1)))
              IF(T1<33.)CYCLE
          ELSE
	          T1=1500.
	      END IF
	      WP2=LOG10(S15(ISL,ISA))
          FC2=LOG10(FC(ISL,ISA))
          T2=10.**(3.1761-1.6576*(LOG10(ST(ISL,ISA))-WP2)/(FC2-WP2))
          IF(T1<T2)CYCLE
          X1=XX*MIN(PRMT(61),(T1-T2)/T1,PKRZ(L1))
          ST(L1,ISA)=ST(L1,ISA)+X1
          PKRZ(L1)=PKRZ(L1)-X1
          ST(ISL,ISA)=MAX(1.E-5,ST(ISL,ISA)-X1)
      END DO
      IF(PKRZ(L1)<0.)PKRZ(L1)=0.
      FPF(ISA)=0.
      RETURN
      END

