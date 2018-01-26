	SUBROUTINE EHYD
!     APEX0806
!     THIS SUB PROGRAM COMPUTES SEDIMENT GRAPHS USING HYDROGRAPH FLOW
!     RATES AND AVE SED CONC
      USE PARM
      QQ=HYDV(IDO)
	  YY=YSD(NDRV,IDO)+STY(IDO)
	  X1=0.
	  !X1=STH(IDO)/CMM
	  R1=.5*X1/(QQ+X1)
	  R2=R1*YY
	  R2=MAX(R2,STY(IDO))
	  YX=YY-R2
	  CIN=.1*YX/QQ
	  T1=0.
      WRITE(KW(26),7)IY,MO,KDA,IDO,ISA,NBSA(ISA),WSA(ISA),QQ,YY,R2
      WRITE(KW(26),9)
      Q1=QHY(1,IDO,IHX(1))
	  SUM=0.
 	  XX=180.*DTHY/WSA(ISA)
      XZ=2.*DTHY 
	  DO I=2,NHY(IDO)
	      Q2=QHY(I,IDO,IHX(1))
	      Q3=QHY(I+1,IDO,IHX(1))
	      C=ABS(Q2+(Q3-Q1)/XZ)
     	  YHY(I,IDO)=XX*(Q2+Q1)*C
 	      WRITE(KW(26),13)I,T1,Q2,C,YHY(I,IDO)
   	      SUM=SUM+YHY(I,IDO)
	      Q1=Q2
	  END DO
      RTO=YX/SUM
	  TOT=0.
	  XX=1.E-6*XX
	  Q1=QHY(1,IDO,IHX(1))
	  DO I=2,NHY(IDO)
	      YHY(I,IDO)=RTO*YHY(I,IDO)
	      T1=T1+DTHY
	      TOT=TOT+YHY(I,IDO)
	      Q2=QHY(I,IDO,IHX(1))
	      X1=XX*(Q1+Q2)
	      IF(X1>0.)THEN
	          C=.1*YHY(I,IDO)/X1
	      ELSE
	          C=0.
	      END IF
   	      IF(KFL(26)>0.)WRITE(KW(26),6)I,T1,Q2,C,YHY(I,IDO)
	      Q1=Q2
	  END DO
	  STY(IDO)=YY-TOT
      WRITE(KW(26),10)YY,R2,TOT,STY(IDO)
	  RETURN
    6 FORMAT(6X,I4,F10.2,F10.4,F10.1,F10.5)
    7 FORMAT(//T10,'COMPUTE SED GRAPH'/T10,'DATE(Y-M-D)=',3I4, 2X,&
     &'IDO= ',I8,2X,'ISA= ',I8,2X,'IDSA= ',I8/T10,'WSA= ',F12.3,' ha',2X&
     &,'HYDV= ',F8.3,' mm',2X,'YI= ',F8.3,' t/ha',2X,'STY= ',F8.3,&
     &' t/ha')
    9 FORMAT(//9X,'#',5X,'Th',10X,'QHYm3/s',5X,'CYppm',7X,'YHYt/s')
   10 FORMAT(10X,'YI=',F8.3,2X,'YIS=',F8.3,2X,'YO=',F8.3,2X,'YOS=',F8.3)
   12 FORMAT(5X,A2,4I4,2I2,5F10.2)
   13 FORMAT(6X,I4,6F12.5)
   15 FORMAT(20X,7F7.2)
	END 	 