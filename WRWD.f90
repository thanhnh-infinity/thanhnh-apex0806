      SUBROUTINE WRWD(I,JNT)
!     APEX0806
!     THIS SUBPROGRAM DETERMINES RAINFALL OCCURENCE AND AMOUNT.
!     AMOUNT IS ESTIMATED BY CALLING SUBPROGRAM WRAIN OR FROM A
!     MODIFIED EXPONENTIAL DISTRIBUTION.
      USE PARM 
      IF(JNT==0)THEN
          RN=1.-AUNIF(IDG(1))
          IF(RN>PRW(LW,IWI,MO)+.001)THEN
              RFV0(I)=0.
              LW=1
              RETURN
          END IF
      END IF
      V4=AUNIF(IDG(3))
      IF(ICDP==0)THEN
          R6=RST(3,IWI,MO)/6.
          ZZ=ADSTN(V3,V4)
          RFV0(I)=WRAIN(R6,ZZ,RST,IWI,MO)*PCF(IWI,MO)
          V3=V4
      ELSE
          RFV0(I)=RST(1,IWI,MO)*(-LOG(V4))**EXPK
      END IF
      LW=2
      RETURN
      END