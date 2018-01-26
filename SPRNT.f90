      SUBROUTINE SPRNT
!     APEX0806
!     THIS SUBPROGRAM PREPARES SOIL DATA TABLE FOR OUTPUT, AND CONVERTS
!     WEIGHTS OF MATERIALS TO CONCENTRATION.
      USE PARM 
      XX=0.
      UW(1)=0.
      UW(2)=0.
      UW(3)=0.
      UW(4)=0.      
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          WT1=WT(ISL,ISA)/1000.
          SOIL(1,ISL,ISA)=WPML(ISL,ISA)/WT1
          SOIL(2,ISL,ISA)=WPMA(ISL,ISA)/WT1
          SOIL(3,ISL,ISA)=WPMS(ISL,ISA)/WT1
          SOIL(4,ISL,ISA)=WPO(ISL,ISA)/WT1
          SOIL(5,ISL,ISA)=WNMN(ISL,ISA)/WT1
          SOIL(6,ISL,ISA)=WON(ISL,ISA)/WT1
          SOIL(7,ISL,ISA)=.1*WOC(ISL,ISA)/WT(ISL,ISA)
          DG=(Z(ISL,ISA)-XX)*1000.
          SSF(ISL,ISA)=S15(ISL,ISA)/DG
          UW(1)=UW(1)+FC(ISL,ISA)
          SOIL(9,ISL,ISA)=FC(ISL,ISA)/DG
          SOIL(8,ISL,ISA)=PO(ISL,ISA)/DG
          SOIL(13,ISL,ISA)=BDD(ISL,ISA)*BD(ISL,ISA)
          UW(2)=UW(2)+ST(ISL,ISA)
          SOIL(12,ISL,ISA)=ST(ISL,ISA)/DG
          UW(3)=UW(3)+S15(ISL,ISA)
          UW(4)=UW(4)+PO(ISL,ISA)
          XX=Z(ISL,ISA)
      END DO
      XX=Z(LID(NBSL(ISA),ISA),ISA)*1000.
      DO I=1,4
          UW(I)=UW(I)/XX
      END DO
      RETURN
      END