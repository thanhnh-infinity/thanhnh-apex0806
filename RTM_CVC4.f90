      SUBROUTINE RTM_CVC4
!     THIS SUBPROGRAM ROUTES HYDROGRAPHS THROUGH A REACH USING THE 
!     MUSKINGUM-CUNGE METHOD WITH VARIABLE C WITH 4 POINT ITERATIVE 
!     SOLUTION
      USE PARM
      DIMENSION QMS(MHP),QMSI(MHP) 
      CMS=RWSA(IDN1)/(DTHY*360.)
      IWH=0
      NIT=0
      IF(HYDV(IDN1)>QTH)IWH=1
      IF(KFL(26)>0.AND.IWH>0)THEN
          WRITE(KW(26),4)IY,MO,KDA,IDN1,IDN2,NBSA(IDN2),IDO,RWSA(IDN1),WSA&
          (IDN2),RFPL(IDN2),HYDV(IDN1),QCAP(IDN2),YSD(NDRV,IDN1)
          WRITE(KW(26),1)
      END IF
      ADD=0.
      TOT=0.
      TM=0.
      QPK=0.
      TPK=0.
      SMI=0.
      NBCX=0
	  NBFX=0
	  QMS=0.
	  QMSI=0.
	  I=1
      DO J=1,MHX
          I=I-1
          DO K=1,NPD
              I=I+1
              QMSI(I)=QHY(K,IDN1,IHX(J))
          END DO
      END DO
      QI1=QHY(1,IDN1,IHX(1))
      QO1=QHY(1,IDO,IHX(1))
      DO L=1,NHY(IDN1)
          IF(QMSI(L)>0.)EXIT
          QMS(L)=0.
          IF(KFL(26)>0)WRITE(KW(26),18)NIT,TM,ADD,ADD,ADD,ADD,ADD,ADD,&
          ADD,ADD,QMSI(L),QMS(L)
          TM=TM+DTHY
      END DO
      ZCH=RCHD(IDN2)
      CBW=RCBW(IDN2)
      XL3=RFPL(IDN2)/3.6
      T1=XL3/.1
      XLS=1000.*RFPL(IDN2)*RFPS(IDN2)
      GS=.01
      RCHX(IDN2)=SQRT(RCHS(IDN2))/RCHN(IDN2)
      SSS=SQRT(RCSS(IDN2)*RCSS(IDN2)+1.)
      AI1=0.
      AO1=0.
      SUM=0.
      I=L
      IIX=I
      ! ROUTING LOOP
      DO 
          I1=I-1
          SMI=SMI+QMSI(I)                                                                         
          QI2=QMSI(I)                                                                         
          QO2=0.
          CALL HQDAV(AI2,CBW,QI2,SSS,ZCH,ZI2,CHW,FPW,IDN2)
          Q3=QI1+QI2+QO1
          GS=MAX(.01,Q3/3.)
          A3=AI1+AI2+AO1
          QIA=MAX(.1,QI1+QI2)
          AIA=AI2+AI1
          ! SOLVE FOR QO2 USING NEWTONS METHOD       
          DO IT=1,10                                                                  
              FU=QO2-GS
              ! IF THE DIFFERENCE IN SUCCESIVE APPROXIMATIONS OF QO2 IS 
              ! ACCEPTABLE SOLUTION IS COMPLETE.  GO TO NEXT TIME STEP.                                                                                   
              IF(ABS(FU/GS)<=1.E-3)EXIT
              IF(IT>1)THEN
                  DF=FU-FU1
                  ! IF THE DIFFERENCE IN SUCCESSIVE APPROXIMATIONS OF THE
                  ! CONVERGENCE FUNCTION IS ACCEPTABLE SOLN COMPLETE. 
                  ! NEXT TIME STEP                                                                                        
                  IF(ABS(DF)<1.E-5)EXIT
                  DG=GS-GS1
                  ! IF THE DIFFERENCE IN SUCCESSIVE GUESSES IS ACCEPTABLE
                  ! SOLN IS COMPLETE
                  IF(ABS(DG)<1.E-5)EXIT                                                                      
                  DFDG=DF/DG
                  XDG=FU/DFDG
                  GS1=GS     
                  X1=ABS(XDG)
                  X2=.5*GS
                  IF(X1>X2)XDG=X2*X1/XDG                                                                    
                  GS=MAX(.01,GS-XDG)                                                                 
              ELSE
                  GS1=GS                                                            
                  GS=GS+.01*FU                                                                   
              END IF
              FU1=FU
              CALL HQDAV(AO2,CBW,GS,SSS,ZCH,ZO2,CHW,FPW,IDN2)
              QQQ=.25*(Q3+GS)
              AAA=.25*(A3+AO2)
              ! CALCULATE M (SLOPE OF Q/A)
              XMC=LOG10((QO1+GS)/QIA)/LOG10((AO1+AO2)/AIA)
              IF(XMC>2.)THEN
                  XMC=2.
              ELSE
                  IF(XMC<1.)XMC=1.
              END IF
              ! CALCULATE CELERITY
              CLTY=XMC*QQQ/AAA
              TT=RFPL(IDN2)/(3.6*CLTY)
              XX=QQQ/((CHW+FPW)*CLTY*XLS)
              X0=.5*(1.-XX)
              X1=DTHY/TT
              ! CALCULATE M-C COEFS
              X2=2.*X0
              X3=1.-X0
              C0=X1+2.*X3
              C1=(X1+X2)/C0
              C2=(X1-X2)/C0
              C3=(2.*X3-X1)/C0
              ! CALCULATE OUTFLOW AT TIME 2.
              QO2=MAX(0.,C1*QI1+C2*QI2+C3*QO1)
          END DO
          IF(IT>20)WRITE(KW(1),'(1X,A)')'MUSKINGUM 4 PT DID NOT CONVERGE'
          NIT=NIT+IT
          AI1=AI2
          AO1=AO2
          QMS(I)=QO2
          IF(I>NHY(IDN1))THEN
              QMSI(I)=.9*QMSI(I1)                                              
              IF(QMS(I1)<=1.)EXIT
          END IF
          IF(KFL(26)>0)WRITE(KW(26),18)IT,TM,XMC,CLTY,TT,X0,C0,C1,C2,C3,&
          QI2,QO2
          TM=TM+DTHY
          SUM=SUM+QO2
          I=I+1
          QI1=QI2
          QO1=QO2
      END DO
      NHY(IDO)=I-1
      RQRB(IDO)=QPK/CMS
      TC(IDO)=TPK
      IF(SUM>0.)THEN
          ! CORRECT OUTFLOW HYD VOL TO AGREE WITH INFLOW
          RTO=SMI/SUM
          TM=0.
          DO I=1,NHY(IDO)
              QMS(I)=QMS(I)*RTO
              TM=TM+DTHY
              IF(QMS(I)>=QPK)THEN
                  QPK=QMS(I)
                  TPK=TM
              END IF
              QAV=.25*(QMSI(I1)+QMSI(I)+QMS(I1)+QMS(I))
	          !IF(QAV>0..AND.IHY==0)CALL RTYNP(QAV,ADD,TOT,I)
	          DF=TOT+TDEG-TDEP-ADD
          END DO
          SUM=0.
          TM=0.
          I=0
          ! PLACE DAY 1 OUTFLOW INTO QHY(IHX(1)) 
          DO K=1,NPD
              I=I+1
              QHY(K,IDO,IHX(1))=QMS(I)
              SUM=SUM+QHY(K,IDO,IHX(1))
              TM=TM+DTHY
          END DO
          ! COMPUTE OUTFLOW VOL
          SUM=(SUM-.5*(QHY(1,IDO,IHX(1))+QHY(NPD,IDO,IHX(1))))/CMS
          HYDV(IDO)=SUM
          ! PLACE DAY 2/MHX OUTFLOW INTO QHY(IHX(J))
          DO J=2,MHX
              I=I-1
              DO K=1,NPD
                  I=I+1
                  QHY(K,IDO,IHX(J))=QMS(I)
              END DO
              IF(I>NHY(IDO))EXIT
          END DO
	      RQRB(IDO)=QPK/CMS
	      TC(IDO)=TPK
	      TCAV(IDO)=TCAV(IDO)+TPK
	      !PRAV(IDO)=PRAV(IDO)+RQRB(IDO)
          !IF(RQRB(IDO)>PRB(IDO))PRB(IDO)=RQRB(IDO)
          !NQRB(IDO)=NQRB(IDO)+1
          IF(TC(IDO)>TCMX(IDO))THEN
              TCMX(IDO)=TC(IDO)
          ELSE
              IF(TC(IDO)<TCMN(IDO))TCMN(IDO)=TC(IDO)
          END IF    
          XO=NHY(IDO)-IIX
          XIT=NIT
          XO=XIT/XO
          IF(NBCX>0)NBCF(IDN2)=NBCF(IDN2)+1
	      IF(NBFX>0)NBFF(IDN2)=NBFF(IDN2)+1
	      NBCT(IDN2)=NBCT(IDN2)+NBCX
	      NBFT(IDN2)=NBFT(IDN2)+NBFX
      ELSE
          NHY(IDO)=0	      
	  END IF
      IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),8)IDN1,IDN2,NBSA(IDN2),SUM,&
      QPK,TPK,XO,TOT,TDEG,TDEP,ADD
      RETURN
    1 FORMAT(//8X,'IT',5X,'Th',8X,'M',9X,'CLTYm/s',3X,'TTh',7X,'X0',8X,&
      'C0',8X,'C1',8X,'C2',8X,'C3',8x,'QIm3/s',4x,'QOm3/s')
    3 FORMAT(T10,'VOL=',F8.3,' mm',3X,'PEAK=',F8.3,' m3/s')
    4 FORMAT(//T10,'ROUTE'/T10,'DATE(Y-M-D)=',3I4,2X,'IDN1= ',I8,2X,&
      'ISA= ',I8,2X,'IDSA= ',I8,2X,'IDO= ',I8/T10,'WSAI= ',F12.3,' ha',&
      2X,'WSA= ',F12.3,' ha',2X,'RFPL= ',F12.3,' km'/T10,'HYD VOL= ',&
      F8.3,' mm',2X,'QCAP= ',F8.3,' m3/s'/T10,'YIN= ',F8.3,' t/ha')    
    6 FORMAT(6X,2I4,20F10.3)
    7 FORMAT(//9X,'#  IT',5X,'Th',8X,'RFmm',6X,'QIm3/s',4X,'Vm/s&
      ',6X,'TTh',7X,'C',9X,'QI+STm3/s',1X,'STm3/s',4X,'QOm3/s',4X,'YIt/s',&
      5X,'DEGt/s',4X'DEPt/s',4X,'YOt/s',5X,'DF')
    8 FORMAT(T10,'IDN1= ',I8,2X,'ISA= ',I8,2X,'IDSA= ',I8,2X,'HYD VOL= ',&
      F7.3,' mm',2X,'PEAK RATE= ',F7.2,' m3/s',2X,'TP= ',F7.2,' h',2X,&
      'AVE IT= ',F6.3/T10,'YI= ',F8.3,' t/ha',2X,'DEG= ',F8.3,' t/ha',&
      2X,'DEP= ',F8.3,' t/ha',2X,'YO= ',F8.3,'t/ha')                                                                                                                                                                                                                        
   17 FORMAT(3I4,F10.3,F10.1,F10.4)
   18 FORMAT(6X,I4,20F10.3)                                                                                                                     
      END      