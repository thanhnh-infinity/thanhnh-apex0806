      SUBROUTINE HRFIN
!     APEX0806
!     THIS SUBPROGRAM READS RAINFALL AT TIME INTERVAL DTHY.
      USE PARM
      DATA IDA1/0/
      I=1 
      RFDT(1)=0.
      DO
	      I=I+1
	   	  READ(KR(30),104,IOSTAT=NFL)IYZ,MOZ,IDZ,THZ,RFDT(I)
		  IF(NFL/=0)EXIT
		  IF(IYZ/=IYR.OR.MOZ/=MO.OR.IDZ/=KDA)EXIT
      END DO
	  NRF=I-1
	  RETURN
  104 FORMAT(3I4,2F10.0)
	  END