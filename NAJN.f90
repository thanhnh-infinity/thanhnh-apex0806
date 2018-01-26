      SUBROUTINE NAJN(UU,AN,DMD,SUX,AJF,JNT)
!     APEX0806
!     THIS SUBPROGRAM COMPUTES ACTUAL N PLANT UPTAKE FROM EACH
!     LAYER(UPTAKE = MINIMUM OF PLANT DEMAND AND SOIL SUPPLY).
      USE PARM
      DIMENSION UU(MSL),AN(MSL,ISA)
      IF(JNT==0)THEN
          SUM=0.
          RT=DMD/(SUX+1.E-20)
          IF(RT>1.)GO TO 4
          DO J=1,LRD(ISA)
              ISL=LID(J,ISA)
              UU(ISL)=UU(ISL)*RT
              SUM=SUM+UU(ISL)
          END DO
          SUX=SUM
          RETURN
      END IF
      SUM=0.
    4 RT=AJF*(DMD-SUX)
      RT1=RT
      DO J=1,LRD(ISA)
          ISL=LID(J,ISA)
          XX=UU(ISL)+RT
          IF(XX<AN(ISL,ISA))GO TO 6
          IF(AN(ISL,ISA)<1.E-10)THEN
              UU(ISL)=0.
              CYCLE
          END IF
          RT=RT-AN(ISL,ISA)+UU(ISL)
          UU(ISL)=AN(ISL,ISA)
          SUM=SUM+UU(ISL)
      END DO
      SUX=SUM
      RETURN
    6 UU(ISL)=XX
      SUX=SUX+RT1
      RETURN
      END