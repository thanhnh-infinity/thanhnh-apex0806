      SUBROUTINE HVOLQ(IVR)
!     APEX0806
!     THIS SUBPROGRAM PREDICTS DAILY RUNOFF VOLUME AND PEAK RUNOFF RATE
!     GIVEN DAILY PRECIPITATION AND SNOW MELT.
      USE PARM
	  IF(IVR==0)THEN
          SUM=0.
          ADD=0.
	      LD1=LID(1,ISA)
          IF(LUN(ISA)==35)THEN
              SCN=25400./CN0(ISA)-254.
              CN=CN0(ISA)
              GO TO 20
          END IF
          IF(NVCN(ISA)==0)THEN
              XX=0.
              DO JJ=1,NBSL(ISA)
                  ISL=LID(JJ,ISA)
                  IF(Z(ISL,ISA)>1.)GO TO 3
                  ZZ=(Z(ISL,ISA)-XX)/Z(ISL,ISA)
                  SUM=SUM+(ST(ISL,ISA)-S15(ISL,ISA))*ZZ/(FC(ISL,ISA)-S15(ISL,ISA))
                  ADD=ADD+ZZ
                  XX=Z(ISL,ISA)
              END DO 
              GO TO 4
            3 ZZ=1.-XX
              SUM=SUM+(ST(ISL,ISA)-S15(ISL,ISA))*ZZ/(FC(ISL,ISA)-S15(ISL,ISA))
              ADD=ADD+ZZ
              GO TO 4
          END IF
          IF(NVCN(ISA)>1)GO TO 19
          DO JJ=1,NBSL(ISA)
              ISL=LID(JJ,ISA)
              IF(Z(ISL,ISA)>1.)GO TO 26
              SUM=SUM+ST(ISL,ISA)-S15(ISL,ISA)
              ADD=ADD+FC(ISL,ISA)-S15(ISL,ISA)
              L1=ISL
          END DO
        !     SCN=MAX(2.,ADD-SUM)**PRMT(37)
          GO TO 4
       26 RTO=(1.-Z(L1,ISA))/(Z(ISL,ISA)-Z(L1,ISA))
          SUM=SUM+(ST(ISL,ISA)-S15(ISL,ISA))*RTO
          ADD=ADD+(FC(ISL,ISA)-S15(ISL,ISA))*RTO
        4 SUM=SUM/ADD
          IF(SUM>0.)GO TO 13
    !     SCN=SMX(ISA)
       31 SCN=SMX(ISA)*(1.-SUM)**2
          GO TO 27
       13 SUM=100.*SUM
          SCN=SMX(ISA)*(1.-SUM/(SUM+EXP(CNSC(1,ISA)-CNSC(2,ISA)*SUM)))
       27 X1=PRMT(15)*(1.-RSD(LD1,ISA))
          IF(X1>-10.)THEN
              SCN=MAX(3.,(SCN-SMX(ISA))*EXP(X1)+SMX(ISA))  
          ELSE
              SCN=SMX(ISA)
          END IF
          IF(STMP(LID(2,ISA),ISA)<-1.)THEN
              SCN=SCN*PRMT(22)
          ELSE
              IF(RFV(IRF(ISA))>0.)SCN=SCN*EXP(PRMT(25)*(.2-AL5))
	      END IF
          CN=25400./(SCN+254.)
          IF(ISCN==0)THEN
              UPLM=MIN(99.5,CN+5.)
              BLM=MAX(1.,CN-5.)
              CN=ATRI(BLM,CN,UPLM,8)
          END IF
          SCN=25400./CN-254.
          IF(SCN>3.)GO TO 20
          SCN=3.
          CN=25400./(SCN+254.)
          GO TO 20
       19 IF(NVCN(ISA)==2)THEN
              DO JJ=1,NBSL(ISA)
                  ISL=LID(JJ,ISA)
                  SUM=SUM+ST(ISL,ISA)-S15(ISL,ISA)
                  ADD=ADD+FC(ISL,ISA)-S15(ISL,ISA)
                  IF(Z(ISL,ISA)>1.)EXIT
              END DO
              SUM=SUM/ADD
              IF(SUM<0.)GO TO 31
              RTO=MIN(.98,SUM)
              SCN=SMX(ISA)*(1.-RTO)
              GO TO 27
          END IF
          IF(NVCN(ISA)==3)THEN
              SCN=25400./CN0(ISA)-254.
              CN=CN0(ISA)
              GO TO 20
          END IF
          SCN=SCI(ISA)
          GO TO 27
       20 SCN2=PRMT(20)*SCN
          TOT=100.
          DO I=1,9
              TOT=TOT-5.
              IF(CN>TOT)EXIT
          END DO
          CNDS(I,ISA)=CNDS(I,ISA)+1.
          IF(RFV(IRF(ISA))<1.E-5)GO TO 11
          SELECT CASE(INFL)
	          CASE(1)
	              X1=RFV(IRF(ISA))-SCN2
                  IF(X1<=0.)GO TO 8
                  QVOL(IDO)=X1*X1/(RFV(IRF(ISA))+(1.-PRMT(20))*SCN)
	          CASE(2,3)
	              CALL HREXP
	          CASE(4)
	              CALL HRUNF
              CASE(5)
                  CALL HRFDTQ
	      END SELECT
	  ELSE
          QVOL(IDO)=EFI(ISA)*RFV(IRF(ISA))
	  END IF
      IF(IFD(ISA)>0)THEN
          IF(DHT(ISA)>.1)THEN
              CALL HFURD
              IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),17)ISA,NBSA(ISA),IYR,&
              MO,KDA,DHT(ISA),RHTT(ISA),QVOL(IDO),DVOL,XHSM(ISA)
              IF(QVOL(IDO)<DVOL-MAX(0.,ST(LD1,ISA)-PO(LD1,ISA)))THEN
                  QVOL(IDO)=0.
                  IF(DHT(ISA)/DKHL(ISA)<.7)DHT(ISA)=DKHL(ISA)
                  GO TO 11
              END IF
          END IF
          DHT(ISA)=0.
          IF(IDRL(ISA)==0.AND.CPHT(JJK,ISA)<1.)DHT(ISA)=DKHL(ISA)
      END IF
    8 IF(URBF(ISA)>0.)THEN
          X1=RFV(IRF(ISA))-2.67
          IF(X1>0.)QURB(IDO)=X1*X1/(RFV(IRF(ISA))+10.69)
          QVOL(IDO)=QVOL(IDO)*(1.-URBF(ISA))+QURB(IDO)*URBF(ISA)
      END IF 
      IF(ITYP>0)THEN
          TC(IDO)=TCC(ISA)+TCS(ISA)/SQRT(RFV(IRF(ISA)))
          QRB=MIN(.95,SCN2/RFV(IRF(ISA)))
          CALL HTR55
      ELSE
          X2=QVOL(IDO)/DUR
          IF(X2>1.)THEN
              X2=X2**.25
          ELSE
              X2=1.
          END IF
          X4=MIN(SPLG(ISA)/360.,TCS(ISA)/X2)
          TC(IDO)=X4+TCC(ISA)/X2
          ALTC=1.-EXP(-TC(IDO)*PRFF)
          X1=ALTC*QVOL(IDO)
          RQRB(IDO)=X1/TC(IDO)
	      QPR(IDO)=X1/X4
      END IF
   11 X1=WSA(ISA)
      IF(DALG(ISA)>0.)QVOL(IDO)=QVOL(IDO)*(X1-DALG(ISA))/X1
      IF(LUN(ISA)==35)ES=RFV(IRF(ISA))-QVOL(IDO)
      RETURN
   17 FORMAT(1X,2I8,1X,3I4,12X,'DHT = ',F6.0,'mm',2X,'RHTT = ',F6.0,&
      'mm',2X,'Q= ',F6.1,'mm',2X,'DV = ',F7.1,'mm',2X,'HUSC = ',F5.2)
      END