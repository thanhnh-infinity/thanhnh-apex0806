      SUBROUTINE SOCIOA(KK)
!     APEX0806
!     THIS SUBPROGRAM OUTPUTS THE SOIL ORGANIC C AND N VARIABLES ANNUALLY
      USE PARM 
      WRITE(KW(21),30)IYR,MO1,KK,ISA,NBSA(ISA)
      WRITE(KW(21),29)TITSO(ISA)
      WRITE(KW(21),29)TITOP(ISA)
      WRITE(KW(21),2)(SID(LORG(LID(J,ISA),ISA)),J=1,NBSL(ISA)),SID(11)
      WRITE(KW(21),3)(Z(LID(I,ISA),ISA),I=1,NBSL(ISA))
      WRITE(KW(21),22)(BDP(LID(I,ISA),ISA),I=1,NBSL(ISA))
      WRITE(KW(21),23)(SAN(LID(I,ISA),ISA),I=1,NBSL(ISA))
      WRITE(KW(21),24)(SIL(LID(I,ISA),ISA),I=1,NBSL(ISA))
      WRITE(KW(21),25)(CLA(LID(I,ISA),ISA),I=1,NBSL(ISA))
      WRITE(KW(21),26)(ROK(LID(I,ISA),ISA),I=1,NBSL(ISA))
      WRITE(KW(21),4)(WLS(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLS(ISA)
      WRITE(KW(21),5)(WLM(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLM(ISA)
      WRITE(KW(21),6)(WLSL(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLSL(ISA)
      WRITE(KW(21),7)(WLSC(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLSC(ISA)
      WRITE(KW(21),8)(WLMC(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLMC(ISA)
      WRITE(KW(21),9)(WLSLC(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLSLC(ISA)
      WRITE(KW(21),10)(WLSLNC(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLSLNC(ISA)
      WRITE(KW(21),11)(WBMC(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZBMC(ISA)
      X1=.001*ZHSC(ISA)
      WRITE(KW(21),12)(WHSC(LID(I,ISA),ISA),I=1,NBSL(ISA)),X1
      X1=.001*ZHPC(ISA)
      WRITE(KW(21),13)(WHPC(LID(I,ISA),ISA),I=1,NBSL(ISA)),X1
      X1=.001*ZOC(ISA)
      WRITE(KW(21),14)(WOC(LID(I,ISA),ISA),I=1,NBSL(ISA)),X1
      WRITE(KW(21),15)(WLSN(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLSN(ISA)
      WRITE(KW(21),16)(WLMN(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZLMN(ISA)
      WRITE(KW(21),17)(WBMN(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZBMN(ISA)
      WRITE(KW(21),18)(WHSN(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZHSN(ISA)
      WRITE(KW(21),19)(WHPN(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZHPN(ISA)
      WRITE(KW(21),20)(WON(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZON(ISA)
      WRITE(KW(21),31)(ECND(LID(I,ISA),ISA),I=1,NBSL(ISA))
      WRITE(KW(21),32)(WSLT(LID(I,ISA),ISA),I=1,NBSL(ISA)),ZSLT(ISA)
      RETURN
    2 FORMAT(T52,'SOIL LAYER NO'/T18,16(4X,A4))
    3 FORMAT(T5,'DEPTH(m)',T20,16F8.2)
    4 FORMAT(T5,'WLS(kg/ha)',T20,16F8.1)
    5 FORMAT(T5,'WLM(kg/ha)',T20,16F8.1)
    6 FORMAT(T5,'WLSL(kg/ha)',T20,16F8.1)
    7 FORMAT(T5,'WLSC(kg/ha)',T20,16F8.1)
    8 FORMAT(T5,'WLMC(kg/ha)',T20,16F8.1)
    9 FORMAT(T5,'WLSLC(kg/ha)',T20,16F8.1)
   10 FORMAT(T5,'WLSLNC(kg/ha)',T20,16F8.1)
   11 FORMAT(T5,'WBMC(kg/ha)',T20,16F8.1)
   12 FORMAT(T5,'WHSC(kg/ha)',T20,16F8.1)
   13 FORMAT(T5,'WHPC(kg/ha)',T20,16F8.1)
   14 FORMAT(T5,'WOC(kg/ha)',T20,16F8.1)
   15 FORMAT(T5,'WLSN(kg/ha)',T20,16F8.1)
   16 FORMAT(T5,'WLMN(kg/ha)',T20,16F8.1)
   17 FORMAT(T5,'WBMN(kg/ha)',T20,16F8.1)
   18 FORMAT(T5,'WHSN(kg/ha)',T20,16F8.1)
   19 FORMAT(T5,'WHPN(kg/ha)',T20,16F8.1)
   20 FORMAT(T5,'WON(kg/ha)',T20,16F8.0)
   22 FORMAT(T5,'BD 33KPA(t/m3)',T20,16F8.2)
   23 FORMAT(T5,'SAND(%)',T20,16F8.1)
   24 FORMAT(T5,'SILT(%)',T20,16F8.1)
   25 FORMAT(T5,'CLAY(%)',T20,16F8.1)
   26 FORMAT(T5,'ROCK(%)',T20,16F8.1)
   29 FORMAT(10X,A12)
   30 FORMAT(//T10,3I4/T10,'SA#= ',I8,1X,'ID= ',I8)
   31 FORMAT(T5,'ECND(mmho/cm)',T20,16F8.2)
   32 FORMAT(T5,'WSLT(kg/ha)',T20,16F8.0)
      END