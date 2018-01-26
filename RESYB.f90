      SUBROUTINE RESYB(YI,YO,YRP,DEP,RSYF,RSYB,ISA,KFL,KW,NBSA,MSO,IPR)
!     APEX0604
!     THIS SUBPROGRAM CHECKS THE RESERVOIR SEDIMENT BALANCE AT THE END
!     OF A SIMULATION.
      DIMENSION NBSA(ISA),KFL(MSO+1),KW(MSO)
      DF=RSYB+YI-YO-YRP-DEP-RSYF
      PER=100.*DF/(DEP+1.E-10)
	  IF(ABS(PER)>1..OR.KFL(1)>0)THEN
	      IF(IPR>0)THEN
	          WRITE(KW(1),2)ISA,NBSA(ISA)
	      ELSE
	          WRITE(KW(1),4)ISA,NBSA(ISA)
	      END IF
	      WRITE(KW(1),3)PER,DF,RSYB,YI,DEP,YO,YRP,RSYF
	  END IF
      RETURN
    2 FORMAT(/T10,'RESERVOIR SED BALANCE',2X,'SA#= ',I8,1X,'ID= ',I8)
    3 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'RSYB=',E13.6,2X,&
     &'YI  =',E13.6,2X,'DEP =',E13.6,2X,'YO  =',E13.6/5X,'YRP =',E13.6,&
     &2X,'RSYF=',E13.6)
    4 FORMAT(/T10,'POND SED BALANCE',2X,'SA#= ',I8,1X,'ID= ',I8)
      END