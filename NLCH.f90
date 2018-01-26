      SUBROUTINE NLCH(SN3I)
!     APEX0806
!     THIS SUBPROGRAM ESTIMATES DAILY SOL N LEACHING BY PERCOLATION AND
!     LATERAL SUBSURFACE FLOW FOR ALL LAYERS EXCEPT THE SURFACE LAYER.
      USE PARM
      WNMN(ISL,ISA)=WNMN(ISL,ISA)+SN3I
      WSLT(ISL,ISA)=WSLT(ISL,ISA)+VSLT(ISA)
      WPML(ISL,ISA)=WPML(ISL,ISA)+VAP(ISA)
      WPMU(ISL,ISA)=WPMU(ISL,ISA)+VPU(ISA)
      SSFN=0.
      SN3I=0.
      VAP(ISA)=0.
      VPU(ISA)=0.
      VSLT(ISA)=0.
      SSST=0.
      QSFN=0.
      IF(PKRZ(ISL)>0.)THEN
	      VH=QSF(ISL,ISA)+CPFH(ISL,ISA)
          V=PKRZ(ISL)+SSF(ISL,ISA)+VH
          VP=V/(PO(ISL,ISA)*PRMT(4))
          IF(VP>5.)THEN
              X1=.99
          ELSE
              X1=1.-EXP(-VP)
          END IF
          IF(WNMN(ISL,ISA)>0.)THEN
              VQN=WNMN(ISL,ISA)*X1
              VQNU=WNMU(ISL,ISA)*X1
              VV=(VQN+VQNU)/(V+1.E-5)
              WNMN(ISL,ISA)=MAX(0.,WNMN(ISL,ISA)-VQN)
              WNMU(ISL,ISA)=MAX(0.,WNMU(ISL,ISA)-VQNU)
              SSO3(ISL)=VV*PKRZ(ISL)
              SSFN=VV*SSF(ISL,ISA)
              QSFN=VV*VH
              SN3I=SSO3(ISL)
          END IF
          IF(WSLT(ISL,ISA)>0.)THEN
              X3=WSLT(ISL,ISA)*X1
              WSLT(ISL,ISA)=MAX(1.E-5,WSLT(ISL,ISA)-X3)
              CS1=X3/(V+1.E-5)
              VSLT(ISA)=CS1*PKRZ(ISL)
              SSST=CS1*(SSF(ISL,ISA)+VH)
          END IF 
          IF(WPML(ISL,ISA)>0..OR.WPMU(ISL,ISA)>0.)THEN
              XX=MIN(.75,PKRZ(ISL)/WT(ISL,ISA))
              VAP(ISA)=XX*WPML(ISL,ISA)
              WPML(ISL,ISA)=WPML(ISL,ISA)-VAP(ISA)
              VPU(ISA)=XX*WPMU(ISL,ISA)
              WPMU(ISL,ISA)=WPMU(ISL,ISA)-VPU(ISA)
          END IF
      END IF
      RETURN
      END