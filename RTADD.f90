      SUBROUTINE RTADD
!     APEX0806
!     THIS SUBPROGRAM ADDS SUBAREA OUTPUTS TO ROUTED OUTPUTS TO
!     DETERMINE TOTAL OUTPUT FROM A ROUTING REACH.
      USE PARM
      IDO=IDOT(ICMD)
      IDN1=IDN1T(ICMD)
      IDN2=IDN2T(ICMD)
      Z1=RWSA(IDN1)
      Z2=RWSA(IDN2)
      Z3=Z1+Z2
      RWSA(IDO)=Z3
      RQRB(IDO)=0.
      QVOL(IDO)=(QVOL(IDN1)*Z1+QVOL(IDN2)*Z2)/Z3
      WYLD(IDO)=(WYLD(IDN1)*Z1+WYLD(IDN2)*Z2)/Z3
      RSSF(IDO)=(RSSF(IDN1)*Z1+RSSF(IDN2)*Z2)/Z3
      QRF(IDO)=(QRF(IDN1)*Z1+QRF(IDN2)*Z2)/Z3
      QDR(IDO)=(QDR(IDN1)*Z1+QDR(IDN2)*Z2)/Z3
      QN(IDO)=(QN(IDN1)*Z1+QN(IDN2)*Z2)/Z3
      QC(IDO)=(QC(IDN1)*Z1+QC(IDN2)*Z2)/Z3
      QP(IDO)=(QP(IDN1)*Z1+QP(IDN2)*Z2)/Z3
      QPU(IDO)=(QPU(IDN1)*Z1+QPU(IDN2)*Z2)/Z3
      YSD(NDRV,IDO)=(YSD(NDRV,IDN1)*Z1+YSD(NDRV,IDN2)*Z2)/Z3
      YC(IDO)=(YC(IDN1)*Z1+YC(IDN2)*Z2)/Z3
      YN(IDO)=(YN(IDN1)*Z1+YN(IDN2)*Z2)/Z3
      YP(IDO)=(YP(IDN1)*Z1+YP(IDN2)*Z2)/Z3
      YMNU(IDO)=(YMNU(IDN1)*Z1+YMNU(IDN2)*Z2)/Z3
      YCOU(IDO)=(YCOU(IDN1)*Z1+YCOU(IDN2)*Z2)/Z3
      YNOU(IDO)=(YNOU(IDN1)*Z1+YNOU(IDN2)*Z2)/Z3
      YPOU(IDO)=(YPOU(IDN1)*Z1+YPOU(IDN2)*Z2)/Z3
      IF(ICDT(ICMD-2)==2)THEN
          II=IDOA(ISA)
          X1=Z1/Z3
          SST(IDO)=SST(II)*X1
          TSFN(IDO)=MAX(0.,TSFN(II)*X1)
      ELSE
          SST(IDO)=(SST(IDN1)*Z1+SST(IDN2)*Z2)/Z3
          TSFN(IDO)=(TSFN(IDN1)*Z1+TSFN(IDN2)*Z2)/Z3
      END IF    
      RSFN(IDO)=(RSFN(IDN1)*Z1+RSFN(IDN2)*Z2)/Z3
      QDRN(IDO)=(QDRN(IDN1)*Z1+QDRN(IDN2)*Z2)/Z3
      QRFN(IDO)=(QRFN(IDN1)*Z1+QRFN(IDN2)*Z2)/Z3
      DO K=1,NDP
          QPST(K,IDO)=(QPST(K,IDN1)*Z1+QPST(K,IDN2)*Z2)/Z3
          TSPS(K,IDO)=(TSPS(K,IDN1)*Z1+TSPS(K,IDN2)*Z2)/Z3
          RSPS(K,IDO)=(RSPS(K,IDN1)*Z1+RSPS(K,IDN2)*Z2)/Z3
          YPST(K,IDO)=(YPST(K,IDN1)*Z1+YPST(K,IDN2)*Z2)/Z3
          VARP(2,K,IDO)=QPST(K,IDO)
          SMMP(2,K,MO,IDO)=SMMP(2,K,MO,IDO)+QPST(K,IDO)
          VARP(4,K,IDO)=TSPS(K,IDO)
          SMMP(4,K,MO,IDO)=SMMP(4,K,MO,IDO)+TSPS(K,IDO)
          VARP(5,K,IDO)=YPST(K,IDO)
          SMMP(5,K,MO,IDO)=SMMP(5,K,MO,IDO)+YPST(K,IDO)
          VARP(11,K,IDO)=RSPS(K,IDO)
          SMMP(11,K,MO,IDO)=SMMP(11,K,MO,IDO)+RSPS(K,IDO)
      END DO
      X1=RSSF(IDO)+SST(IDO)+QRF(IDO)+QDR(IDO)
	  WYLD(IDO)=QVOL(IDO)+X1
      IF(WYLD(IDO)>0.)THEN
		  IF(QVOL(IDO)>0.)THEN
              TC(IDO)=(TC(IDN1)*QVOL(IDN1)+TC(IDN2)*QVOL(IDN2))/&
              (QVOL(IDN1)+QVOL(IDN2))
              IF(TC(IDO)>0.)THEN
	              TCAV(IDO)=TCAV(IDO)+TC(IDO)
			      ALTC=1.-EXP(-TC(IDO)*PRFF)
			      RQRB(IDO)=ALTC*QVOL(IDO)/TC(IDO)
			      PRAV(IDO)=PRAV(IDO)+RQRB(IDO)
			      IF(RQRB(IDO)>PRB(IDO))PRB(IDO)=RQRB(IDO)
			      NQRB(IDO)=NQRB(IDO)+1
			  END IF
	      ELSE
		      TC(IDO)=MAX(TC(IDN1),TC(IDN2))
			  PRFF=.042
			  ALTC=1.-EXP(-TC(IDO)*PRFF) 
	      END IF
	      !IF(KFL(9)>0.AND.ABS(Z3-RWSA(NCMD))<1.E-5)THEN
			  !X1=MAX(RQRB(IDO),X1/24.)
              !WRITE(KW(9),1202)IY,IYR,MO,KDA,QVOL(IDO),SST(IDO),&
              !QRF(IDO),RSSF(IDO),WYLD(IDO),X1,TC(IDO),ALTC
		  !END IF
	  END IF
      X1=YSD(NDRV,IDN1)*Z1
      X2=YSD(NDRV,IDN2)*Z2
      X3=X1+X2
      TOT=0.
      DO I=1,NSZ
          PCT(I,IDO)=MAX(.01,(PCT(I,IDN1)*X1+PCT(I,IDN2)*X2)/(X3+1.E-5))
          TOT=TOT+PCT(I,IDO)
      END DO
      PSZM(IDO)=0.
      DO I=1,NSZ
          PCT(I,IDO)=PCT(I,IDO)/(TOT+1.E-10)
          PSZM(IDO)=PSZM(IDO)+PSZY(I)*PCT(I,IDO)
      END DO
      IF(KFL(NOFL)>0.AND.ICMD==ICMO(IOF))THEN
	      IF(IY==1.AND.IDA==1)WRITE(KW(NOFL),4082)Z3
          X2=10.*Z3
          XTP(1)=X2*(QVOL(IDO)+RSSF(IDO)+QRF(IDO)+SST(IDO)+QDR(IDO))
	      SMSO(IOF)=SMSO(IOF)+XTP(1)
          XTP(2)=Z3*YSD(NDRV,IDO)
          XTP(3)=Z3*YN(IDO)
          XTP(4)=Z3*YP(IDO)
          X1=QN(IDO)+RSFN(IDO)+QRFN(IDO)+TSFN(IDO)+QDRN(IDO)
          XTP(5)=Z3*X1
          XTP(6)=Z3*(QP(IDO)+QPU(IDO))
          WRITE(KW(NOFL),902)IDA,IYR,(XTP(I),I=1,6),(PSZ(I),PCT(I,IDO),&
          I=1,NSZ)
	      IOF=IOF+1
	      NOFL=NOFL+1
	  END IF
      IF(IHY==0)RETURN
      NHY(IDO)=MAX(NHY(IDN1),NHY(IDN2))
      IF(NHY(IDO)==0)RETURN
      IWH=0
      IF(QVOL(IDO)>QTH)IWH=1
      IF(KFL(26)>0.AND.IWH>0)THEN
          WRITE(KW(26),11)IY,MO,KDA,IDN1,IDN2,IDO
          WRITE(KW(26),9)
      END IF
      T1=0.
      QPK=0.
      SUM=0.
      DO K=1,NPD
          QHY(K,IDO,IHX(1))=QHY(K,IDN1,IHX(1))+QHY(K,IDN2,IHX(1))
          SUM=SUM+QHY(K,IDO,IHX(1))
          IF(QHY(K,IDO,IHX(1))>QPK)THEN
              QPK=QHY(K,IDO,IHX(1))
              TPK=T1
          END IF
          IF(KFL(26)>0)WRITE(KW(26),6)T1,QHY(K,IDN1,IHX(1)),QHY(K,IDN2,&
          IHX(1)),QHY(K,IDO,IHX(1))
          T1=T1+DTHY
      END DO
      X2=Z3/(DTHY*360.)
      RQRB(IDO)=QPK/X2
      SUM=(SUM-.5*(QHY(1,IDO,IHX(1))+QHY(NPD,IDO,IHX(1))))/X2
      HYDV(IDO)=SUM
      TC(IDO)=TPK
      !NQRB(IDO)=NQRB(IDO)+1
      !PRAV(IDO)=PRAV(IDO)+RQRB(IDO)
      !PRSD(IDO)=PRSD(IDO)+RQRB(IDO)*RQRB(IDO)
      !IF(RQRB(IDO)>PRB(IDO))PRB(IDO)=RQRB(IDO)
      !X2=RQRB(IDO)/(SUM+1.E-10)
      !IF(X2>QRQB(IDO))QRQB(IDO)=X2
      !QRBQ(IDO)=QRBQ(IDO)+X2
      TCAV(IDO)=TCAV(IDO)+TC(IDO)
      IF(TC(IDO)<TCMX(IDO))THEN
          IF(TC(IDO)<TCMN(IDO))TCMN(IDO)=TC(IDO)
      ELSE
          TCMX(IDO)=TC(IDO)
      END IF
      IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),13)Z1,Z2,Z3,HYDV(IDN1),&
      HYDV(IDN2),SUM,QPK,TPK
      I=24
      DO J=2,MHX
          T1=T1-DTHY
          DO K=1,NPD
              I=I+1
              QHY(K,IDO,IHX(J))=QHY(K,IDN1,IHX(J))+QHY(K,IDN2,IHX(J))
              IF(KFL(26)>0)WRITE(KW(26),6)T1,QHY(K,IDN1,IHX(J)),QHY(K,IDN2,&
              IHX(J)),QHY(K,IDO,IHX(J))
              T1=T1+DTHY
          END DO
          IF(I>NHY(IDO).AND.QHY(NPD,IDO,IHX(J))<.1)EXIT
	  END DO
      RETURN
    6 FORMAT(9X,F8.3,4F10.3)
    9 FORMAT(//13X,'Th',5X,'QHY1m3/s',2X,'QHY2m3/s',2X,'QHYOm3/s')
   11 FORMAT(//T10,'ADD HYD'/T10,'DATE(Y-M-D)=',3I4,2X,'IDN1= ',I8,2X&
      ,'IDN2= ',I8,2X,'IDO= ',I8)
   13 FORMAT(T10,'WSA(IDN1)= ',F12.3,' ha',2X,'WSA(IDN2)= ',F12.3,' ha',&
      2X,'WSA(IDO)= ',F12.3,' ha'/T10,'QID1= ',F8.3,' mm',2X,'QID2= ',&
      F8.3,' mm',2X, 'QO= ',F8.3,' mm'/T10,'PEAK RATE= ',F8.3,' m3/s',2X,&
      'TP= ',F7.2,' h')
   14 FORMAT(T10,'ADD SED= ',F7.3,' t/ha')  
   12 FORMAT(5X,A2,3I8,3I4,5F10.2)
   17 FORMAT(20X,10E13.5)
  902 FORMAT(1X,I4,1X,I4,1X,20(1X,E16.10))
 1202 FORMAT(1X,4I4,15F10.2)  
 4082 FORMAT(///10X,'WATERSHED AREA = ',F10.2,' HA'/1X,'JDA   YR     TOT &
      FLOWm3        YSDt             YNkg             YPkg          TOT &
      NO3kg          QPkg',11X,3('PSZum',11X,'FRACTION',8X))  
      END