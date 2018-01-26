      SUBROUTINE AISHFL
!     APEX0806
!     THIS SUBPROGRAM SHUFFLES DATA RANDOMLY.(BRATLEY,FOX,SCHRAGE,P.34)
      USE PARM
      DO I=12,2,-1
          II=IDG(I)
          RN=AUNIF(13)
          K=I*RN+1
          IDG(I)=IDG(K)
          IDG(K)=II
      END DO
      RETURN
      END