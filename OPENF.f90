      SUBROUTINE OPENF(ASTN)
!     APEX0806
!     THIS SUBPROGRAM OPENS FILES.
      USE PARM
      CHARACTER(4)::AXT
      CHARACTER(20)::ASTN
      DIMENSION AXT(44)
	  DATA AXT/".OUT",".MAN",".SUS",".ASA",".SWT",".DPS",".MSA",".AWP",&
      ".DHY",".WSS",".SAD",".HYC",".DRS","    ",".MWS",".DWS",".AWS",&
      ".DGZ",".DUX",".DDD",".ACN",".DCN",".SCX",".ACY",".EFR",".EHY",&
      ".APS",".MSW",".DPW",".SPS",".ACO",".SWN",".FSA",".SAO",".RCH",&
      ".ERX",".MRH",".STR","    ",".SW1",".SW2",".SW3",".SW4",".SW5"/
      OPEN(KW(1),FILE=ASTN//AXT(1))
	  DO I=2,44
	      IF(AXT(I)/="    ".AND.KFL(I)>0)OPEN(KW(I),FILE=ASTN//AXT(I))
	  END DO
      RETURN
      END
      