      SUBROUTINE OPENV(NUM,FNAM,IDIR)
!     VERIFIES THE EXISTENCE OF A FILE BEFORE OPENING IT
      CHARACTER(20)::FNAM,DIR
	  CHARACTER(36)::FNM
	  LOGICAL::XMIS	
	  DIR='C:\WEATDATA\'
	  IF(IDIR/=0)THEN
          FNM=ADJUSTR(DIR)//ADJUSTL(FNAM)
	  ELSE
	      FNM=FNAM
	  END IF	
	  INQUIRE(FILE=FNM,EXIST=XMIS)
	  IF(XMIS==.TRUE.)THEN
	      OPEN(NUM,FILE=FNM)
	  ELSE
          WRITE(*,'(/A/)')'File '//TRIM(FNM)//' IS MISSING.'
          PAUSE
          STOP
	  END IF	
	  RETURN
      END