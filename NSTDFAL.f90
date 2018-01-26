      SUBROUTINE NSTDFAL(ZZ)
!     APEX0806
!     THIS SUBPROGRAM SIMULATES THE CONVERSION OF STANDING DEAD CROP
!     RESIDUE TO FLAT RESIDUE.
      USE PARM
      LD1=LID(1,ISA)
      SUM=0.
      TOT=0.
      DO K=1,LC
          IF(STD(K,ISA)<.001)CYCLE
          Z1=ZZ*STD(K,ISA)
          STD(K,ISA)=STD(K,ISA)-Z1
          SUM=SUM+Z1
          Z2=ZZ*STDN(K,ISA)
          STDN(K,ISA)=MAX(1.E-5,STDN(K,ISA)-Z2)
          TOT=TOT+Z2
          Z3=ZZ*STDP(K,ISA)
          STDP(K,ISA)=STDP(K,ISA)-Z3
          FOP(LD1,ISA)=FOP(LD1,ISA)+Z3
          Z4=ZZ*STDL(K,ISA)
          STDL(K,ISA)=MAX(1.E-5,STDL(K,ISA)-Z4)
      END DO
      ZZ=MIN(1.,ZZ*10.)
      ZS=ZZ*STDO(ISA)
      STDO(ISA)=MAX(1.E-5,STDO(ISA)-ZS)
      SUM=SUM+ZS
      ZS=ZZ*STDON(ISA)
      TOT=TOT+ZS
      CALL NCNSTD(.05,SUM,TOT,LD1)
      STDON(ISA)=MAX(1.E-5,STDON(ISA)-ZS)
      ZS=ZZ*STDOP(ISA)
      FOP(LD1,ISA)=FOP(LD1,ISA)+ZS
      STDOP(ISA)=MAX(1.E-5,STDOP(ISA)-ZS)
      RETURN
      END