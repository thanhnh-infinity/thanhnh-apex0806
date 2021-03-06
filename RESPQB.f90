      SUBROUTINE RESPQB(QI,EV,QO,SP,QRP,VRSE,RSVB,ISA,KFL,KW,NBSA,MSO,IPR)
!     APEX0604
!     THIS SUBPROGRAM CHECKS THE POND AND RESERVOIR WATER BALANCES AT 
!     THE END OF A SIMULATION
      DIMENSION NBSA(ISA),KFL(MSO+1),KW(MSO)
      DF=RSVB+QI-EV-QO-VRSE-SP-QRP
	  PER=200.*DF/(QI+QO)
	  IF(ABS(PER)>1..OR.KFL(1)>0)THEN
	      IF(IPR>0)THEN
	          WRITE(KW(1),4)ISA,NBSA(ISA)
	      ELSE
	          WRITE(KW(1),2)ISA,NBSA(ISA)
	      END IF
	      WRITE(KW(1),3)PER,DF,RSVB,QI,QO,EV,SP,QRP,VRSE
          RSVB=VRSE
      END IF
      RETURN
    2 FORMAT(/T10,'POND WATER BALANCE',2X,'SA#= ',I8,1X,'ID= ',I8)
    3 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'BSW =',E13.6,2X,&
      'QI  =',E13.6,2X,'QO  =',E13.6/5X,'EVP =',E13.6,2X,'SEP =',E13.6,&
      2X,'QRP =',E13.6,2X,'FSW =',E13.6)
    4 FORMAT(/T10,'RESERVOIR WATER BALANCE',2X,'SA#= ',I8,1X,'ID= ',I8)
      END