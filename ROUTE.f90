      SUBROUTINE ROUTE
!     APEX0806
!     THIS SUBPROGRAM CONTROLS ROUTING OPERATIONS.
      USE PARM 
      IDO=IDOT(ICMD)
      IDN1=IDN1T(ICMD)
      IDN2=IDN2T(ICMD)
      IDNB(IDO)=NBSA(IDN2)
      QVOL(IDO)=QVOL(IDN1)
      RSSF(IDO)=RSSF(IDN1)
      SST(IDO)=SST(IDN1)
      QDR(IDO)=QDR(IDN1)
      QRF(IDO)=QRF(IDN1)
      WYLD(IDO)=QVOL(IDO)+RSSF(IDO)+SST(IDO)+QDR(IDO)+QRF(IDO)
      RWSA(IDO)=RWSA(IDN1)
      TC(IDO)=TC(IDN1)
      RSFN(IDO)=RSFN(IDN1)
      QDRN(IDO)=QDRN(IDN1)
      QRFN(IDO)=QRFN(IDN1)
      RQRB(IDO)=RQRB(IDN1)
      DO I=1,NSZ
          PCT(I,IDO)=PCT(I,IDN1)
      END DO
      DO K=1,NDP
          QPST(K,IDO)=QPST(K,IDN1)
          TSPS(K,IDO)=TSPS(K,IDN1)
	      YPST(K,IDO)=YPST(K,IDN1)
      END DO
      X2=RWSA(IDN1)/WSA(IDN2)
      X6=SST(IDN1)*X2
      X5=TSFN(IDN1)*X2
      SST(IDO)=SST(IDN1)
      II=IDOA(IDN2)
      SST(II)=SST(IDN1)
      TSFN(IDO)=TSFN(IDN1)
      TSFN(II)=TSFN(IDN1)
      SSIN(IDN2)=SSIN(IDN2)+X5
      X2=X5/Z(LID(NBSL(IDN2),IDN2),IDN2)
      SSFI(IDN2)=SSFI(IDN2)+X6
      SMM(95,MO,IDN2)=SMM(95,MO,IDN2)+X6
      VAR(95,IDN2)=X6
      XX=X6/Z(LID(NBSL(IDN2),IDN2),IDN2)
      Z1=0.
      DO I=1,NBSL(IDN2)
          ISL=LID(I,IDN2)
          X1=Z(ISL,IDN2)-Z1
          X3=XX*X1
          X4=X2*X1
          WNMN(ISL,IDN2)=WNMN(ISL,IDN2)+X4
          ST(ISL,IDN2)=ST(ISL,IDN2)+X3
          Z1=Z(ISL,IDN2)
      END DO
      YSD(NDRV,IDO)=0.
      YN(IDO)=0.
      YP(IDO)=0.
      YC(IDO)=0.
      QC(IDO)=0.
      QN(IDO)=0.
      QP(IDO)=0.
	  QPU(IDO)=0.
      YMNU(IDO)=0.
      YCOU(IDO)=0.
      YNOU(IDO)=0.
      YPOU(IDO)=0.
      IF(IHY==0.AND.WYLD(IDO)>0.)THEN
          CALL RTSED
      ELSE    
          IF(NHY(IDN1)>0.)THEN
              SELECT CASE(IHY)
                  CASE(1)
                      CALL RTVSC
                  CASE(2)
                      CALL RTSVS
                  CASE(3)
                      CALL RTM_CVC
                  CASE(4)
                      CALL RTM_CVC4
                  CASE DEFAULT
                      CALL RTVSC
              END SELECT                                      
          END IF
      END IF    
      RETURN
      END