      SUBROUTINE SWNN(CL,SA,OC,W1,F3)
!     APEX0806
!     THIS SUBPROGRAM ESTIMATES FC AND WP USING NEAREST NEIGHBOR 
!     APPROACH.
      USE PARM
      DIMENSION DXS(NSX),XTX(3)
      XTX(1)=(SA-XAV(1))/XDV(1)
      XTX(2)=(CL-XAV(2))/XDV(2)
      XTX(3)=(OC-XAV(3))/XDV(3)
      XTX(1)=XTX(1)*BRNG/XRG(1)
      XTX(2)=XTX(2)*BRNG/XRG(2)
      XTX(3)=XTX(3)*BRNG/XRG(3)
      I1=1
      DO I=1,NSX
          DXS(I)=0.
          DO K=1,3
              DXS(I)=DXS(I)+(XSP(I,K)-XTX(K))**2
          END DO
          DXS(I)=SQRT(DXS(I))
          KX(I1)=I
          I1=I1+1
      END DO
      I1=I1-1
      CALL ASORT3(DXS,KX,NSX)
      N1=MIN(I1,NSNN)
      SUM=0.
      DO I=1,N1
          X1=.01*XSP(KX(I),4)
          X3=X1+XSP(KX(I),5)
          DO K=1,3
              XTX(K)=XDV(K)*XSP(KX(I),K)*XRG(K)/BRNG+XAV(K)
          END DO
          SUM=SUM+DXS(I)
      END DO
      TOT=0.
      DO I=1,N1
          DXS(I)=(SUM/DXS(I))**EXNN
          TOT=TOT+DXS(I)
      END DO
      W1=0.
      F3=0.
      DO I=1,N1
          DXS(I)=DXS(I)/TOT
          X1=.01*XSP(KX(I),4)
          X3=X1+XSP(KX(I),5)
          W1=W1+DXS(I)*X1
          F3=F3+DXS(I)*X3
      END DO
      RETURN
    1 FORMAT(10F10.4)
      END