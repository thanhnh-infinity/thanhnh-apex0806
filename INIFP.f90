      SUBROUTINE INIFP(IDOF,I3,II,JJ,JRT)
!     APEX0806
!     THIS SUBPROGRAM ALLOWS INPUT TO OPERATION SCHEDULE FOR IRRIGATION,
!     FERTILIZER, OR PESTICIDE APPLICATION
      USE PARM
      DIMENSION IDOF(MNT)
      I1=I3-6
      SELECT CASE(I1)
          CASE(1)
              KP=KP+1
              NPC(KP)=JJ
              KPC(KP)=JX(2)
              JPC(KP)=JX(3)
              PSTR(II,KP,ISA)=OPV(1)
              PSTE(II,KP,ISA)=OPV(2)
              CALL PSTTBL
              LPC(II,KP,ISA)=KDP1(JX(7))
          CASE(2)
              VIRR(II,JJ,ISA)=OPV(1)
              IF(ABS(OPV(3))>1.E-5)BIR(ISA)=OPV(3)
              IF(OPV(4)>0.)EFI(ISA)=OPV(4)
              IF(OPV(5)>0.)FIRG(ISA)=OPV(5)
              KI=KI+1
              NIR(KI)=JX(2)
              IIR(KI)=JX(3)
              KIR(KI)=JJ
          CASE(3)
              KF=KF+1
              IHT(KF,ISA)=JJ
              IDOF(KF)=JJ
              NFT(KF)=JX(2)
              IFT(KF)=JX(3)
              WFA(II,KF,ISA)=OPV(1)
              FDP(II,KF,ISA)=OPV(2)*.001
              CALL NFTBL(L)
              LFT(II,KF,ISA)=KDF1(JX(7))      
          CASE DEFAULT
              GO TO 169
      END SELECT
      TIR(II,JJ,ISA)=BIR(ISA)
      CND(II,JJ,ISA)=CN2(ISA)
      QIR(II,JJ,ISA)=EFI(ISA)
      FIRX(II,JJ,ISA)=FIRG(ISA)
      JRT=1
      RETURN
  169 IF(I3==NHC(19))THEN
          RSTK(II,JJ,ISA)=OPV(1)
          IF(OPV(1)>0.)IGZX(1,1)=ISA
      END IF
      IF(OPV(2)<0.)THEN
          CN2(ISA)=-OPV(2)      
      ELSE
          IF(OPV(2)>0.)THEN
	          LUN(ISA)=OPV(2)
	          LUN(ISA)=LUN(ISA)+LUNS(ISA)
              CALL HSGCN
          END IF
      END IF
      CND(II,JJ,ISA)=CN2(ISA)
      IF(ABS(OPV(3))>1.E-5)BIR(ISA)=OPV(3)
      TIR(II,JJ,ISA)=BIR(ISA)
      IF(OPV(4)>0.)EFI(ISA)=OPV(4)
      QIR(II,JJ,ISA)=EFI(ISA)
      FIRX(II,JJ,ISA)=FIRG(ISA)
      JRT=0
      RETURN
      END