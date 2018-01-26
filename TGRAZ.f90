      SUBROUTINE TGRAZ(JRT)
!     APEX0806
!     THIS SUBPROGRAM SIMULATES ANIMAL GRAZING.
      USE PARM 
	  JRT=0
	  LGZ=1
      II=IHC(JT1)
      LD1=LID(1,ISA)
      NN=NCP(IRO(ISA),ISA)
      NGZ=0
	  YY=0.
      YLN=0.
      YLP=0.
      IOW=IDON(ISA)
	  KOMP(KT(ISA),ISA)=0
	  X1=RSTK(IRO(ISA),KT(ISA),ISA)
	  IF(X1>0.)THEN
	      IHD=IHDM(ISA)
	      GCOW(IHD,ISA)=WSA(ISA)/X1
	      IGZ(ISA)=1
	  ELSE
	      DO IHD=1,NHRD(IOW)
              IF(IGZX(IHD,IOW)==ISA)GO TO 2
          END DO
          RETURN
      END IF          
    2 NGZ(IHD,ISA)=IHD
      GZLX=GZLM(IHD,ISA)
	  IF(ORHI(JT1)>0.)GZLX=ORHI(JT1)
      IF(AGPM(ISA)<GZLX.OR.GCOW(IHD,ISA)<1.E-5)THEN
          GCOW(IHD,ISA)=0.
          NGZ(IHD,ISA)=0
          RETURN
	  END IF
      DO K=1,NN
          JJ=JE(K,ISA)
          IF(IDC(JJ)==NDC(7).OR.IDC(JJ)==NDC(8).OR.IDC(JJ)==NDC(10))CYCLE
          IF(JJ==12)CYCLE
          IF(JP(JJ,ISA)<=0)THEN
              JP(JJ,ISA)=1
              NCR(JJ,ISA)=NCR(JJ,ISA)+1
          END IF
          HUF(JJ,ISA)=MAX(HUF(JJ,ISA),HU(JJ,ISA))
          DMF(JJ,ISA)=DM1(JJ,ISA)
          TRA(JJ,ISA)=SRA(JJ,ISA)+TRA(JJ,ISA)
          IF(RD(JJ,ISA)>RDF(JJ,ISA))RDF(JJ,ISA)=RD(JJ,ISA)
          XX=DM(JJ,ISA)+.001
	      X3=RW(JJ,ISA)
	      X12=STL(JJ,ISA)
          X2=UN1(JJ,ISA)/XX
          X7=X2*.001
          X10=UP1(JJ,ISA)/XX
          X8=STD(JJ,ISA)+STDO(ISA)+1.E-10
          RNR=STDN(JJ,ISA)/X8
          RPR=STDP(JJ,ISA)/X8
	      XX=GCOW(IHD,ISA)*GZRT(IHD,IOW)/WSA(ISA)
          X1=MIN(XX/AGPM(ISA),.9)
          AJHI(JJ,ISA)=0.
	      ZZ=MAX(.01,1.-X1)
          CPHT(JJ,ISA)=MAX(.001,CPHT(JJ,ISA)*ZZ)
          HU(JJ,ISA)=MAX(.1,HU(JJ,ISA)*ZZ)
          SLAI(JJ,ISA)=MAX(.01,SLAI(JJ,ISA)*ZZ)
          STD(JJ,ISA)=MAX(.01,STD(JJ,ISA)*ZZ)
          STDO(ISA)=MAX(.01,STDO(ISA)*ZZ)
          YZ=X1*X8      
          XZ=X1*X12  
          YLD(JJ)=XZ*HE(JT1)
          YLSD=YZ*HE(JT1)
          Y4=YZ*RNR
          Y5=YZ*RPR
          STDN(JJ,ISA)=MAX(1.E-10,STDN(JJ,ISA)-Y4)
          STDP(JJ,ISA)=MAX(1.E-10,STDP(JJ,ISA)-Y5)
          X4=MIN(XZ*X2,UN1(JJ,ISA))
          X5=MIN(XZ*X3,UP1(JJ,ISA))
          Z2=YLSD*RNR
          Z3=YLSD*RPR
          YLN=MIN(.9*(UN1(JJ,ISA)+STDN(JJ,ISA)),YLD(JJ)*X2)
          YLP=MIN(.9*(UP1(JJ,ISA)+STDP(JJ,ISA)),YLD(JJ)*X3)
          X11=XZ-YLD(JJ)+YZ-YLSD
          X10=X4-YLN+Y4-Z2
	      JJK=JJ
          CALL NCNSTD(.05,X11,X10,LD1)
          FOP(LD1,ISA)=MAX(.01,FOP(LD1,ISA)+X5-YLP+Y5-Z3)
          YY=YLD(JJ)+YLSD
          YLD2(JJ,ISA)=YLD2(JJ,ISA)+YY
          JD(ISA)=JJ
          SRA(JJ,ISA)=0.
          UN1(JJ,ISA)=UN1(JJ,ISA)-X4
          UP1(JJ,ISA)=UP1(JJ,ISA)-X5
          DM(JJ,ISA)=DM(JJ,ISA)-XZ
	      IF(DM(JJ,ISA)<RW(JJ,ISA))RW(JJ,ISA)=DM(JJ,ISA)
          STL(JJ,ISA)=DM(JJ,ISA)-RW(JJ,ISA)
          YLN=YLN+Z2
          YLP=YLP+Z3
          YLNF(JJ,ISA)=YLNF(JJ,ISA)+YLN
          YLPF(JJ,ISA)=YLPF(JJ,ISA)+YLP
          TYN(ISA)=TYN(ISA)+YLN
          TYP(ISA)=TYP(ISA)+YLP
          IF(YLD(JJ)>1.E-10)THEN
	          IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),29)ISA,NBSA(ISA),IYR,&
              MO,KDA,IOW,IHD,TIL(JT1),CPNM(JD(ISA)),YLD(JJ),YLSD,AGPM(ISA),X12,&
              X3,X1,X7,XHSM(ISA),YLN,YLP
              IF(KFL(18)>0)WRITE(KW(18),3)ISA,NBSA(ISA),IYR,MO,KDA,IY,IOW,IHD&
              ,TIL(JT1),CPNM(JD(ISA)),YLD(JJ),YLSD,AGPM(ISA),X12,X8,RSD(LD1,ISA)
          END IF 
      END DO          
      RETURN
    3 FORMAT(1X,2I8,1X,I4,2I2,3I4,1X,A8,1X,A4,6F8.4)
   29 FORMAT(1X,2I8,1X,I4,2I2,2X,'IDON=',I4,2X,'HRD#=',I3,2X,A8,2X,A4,2X&
      ,'YLD=',F7.4,'t/ha',2X,'YSD=',F7.4,'t/ha',2X,'AGPM=',F7.4,'t/ha',&
      2X,'STL=',F7.4,'t/ha',2X,'RWT=',F7.4,'t/ha',2X,'HI=',F7.4,2X,&
      'NCN=',F6.3,'G/G',2X,'HUSC=',F5.2,2X,'YLN=',F5.0,'kg/ha',2X,'YLP='&
      ,F5.0,'kg/ha')
      END