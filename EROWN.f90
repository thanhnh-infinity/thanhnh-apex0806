      FUNCTION EROWN(Y1)
!     APEX0806
!     THIS SUBPROGRAM COMPUTES POTENTIAL WIND EROSION RATE FOR BARE
!     SOIL GIVEN WIND SPEED.
      USE PARM
      DU10=U10(IRF(ISA))*Y1
      USTR=.0408*DU10
      X1=USTR*USTR-USTW
      IF(X1>0.)THEN
          YWR=.255*X1**1.5
          EROWN=YWR*ALG
      ELSE
    !     WRITE(6,15)X,DU10,USTR,YWR,WDST(I)
          EROWN=0.
      END IF
      RETURN
      END