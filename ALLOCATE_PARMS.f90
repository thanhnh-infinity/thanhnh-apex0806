      SUBROUTINE ALLOCATE_PARMS                                                      
!     THIS SUBROUTINE ALLOCATES ARRAY SIZES                                          
      USE PARM                                                                       
      MSL=10                                                                              
      MSO=49                                                                         
	  NSM=140                                                                             
	  NSH=35 
	  ML1=MSL+1                                                                           
      MHY=MSA*4                                                                      
      MHY5=MHY*5                                                                     
      MLA=MSL*MSA                                                                    
      MLA1=ML1*MSA                                                                   
      MLA17=MSL*MSA*17                                                               
	  MPA=MPS*MSA                                                                         
      MPA90=MPA*90                                                                   
      MPLA=MPA*MSL                                                                   
	  MCA=MNC*MSA                                                                         
      MCA5=MCA*5                                                                     
      MCA12=MCA*12                                                                   
      MCCA=MNC*MNC*MSA                                                               
      MRA=MRO*MSA                                                                    
      MRCA=MRA*MNC                                                                   
      MRTA=MNT*MRA                                                                   
	  MLCA=MNC*MLA
	  NBMX=MAX(10000,MHY)                                                                        
	  ALLOCATE (CPNM(MNC),FPSO(NBMX),FTNM(MFT),PSTN(MPS),TITOP(MSA),TITSA&
	  (MSA),TITSO(MSA),HEDH(NSH),HED(NSM))                                                            
	  ALLOCATE (KW(2*MSA+MSO),ICDT(MHY),IDN1T(MHY),IDN2T(MHY),IDNB(MHY),&                 
      IDOT(MHY),IDRO(MHY),NHY(MHY),NQRB(MHY),NTX(MHY),NISA(NBMX),ICUS&
      (MNT),IHC(MNT),NBE(MNT),NBT(MNT),IDC(MNC),KDC(MNC),NTP(MNC),IFT&
      (MFT),KDF(MFT),NFT(MFT),KFL(MSO+1))                                                     
	  ALLOCATE (ISAL(MOW),ISAS(MOW),NFED(MOW),NHRD(MOW),NSAL(MSA),NSAO&                   
      (MOW),NSAS(MOW))                                                               
	  ALLOCATE (IIR(MRO),KIR(MRO),NIR(MRO))                                               
	  ALLOCATE (IAC(MSA),IAMF(MSA),IAPL(MSA),IAUF(MSA),IAUI(MSA),IAUL&                    
      (MSA),IBSA(MSA),IDFH(MSA),IDGH(MSA),IDOA(MSA),IDR(MSA),IDRL(MSA),&             
      IDON(MSA),IDS(MSA),IEXT(MSA),IFA(MSA),IFD(MSA),IFLS(MSA),IGO(MSA),&            
      IGZ(MSA),IHDM(MSA),ILQF(MSA),IMW(MSA),IPMP(MSA),IPSO(MSA),IPST(MSA),&
      IPTS(MSA),IRF(MSA),IRI(MSA),IRO(MSA),IRR(MSA),ISAO(MSA),ISCP(MSA))
      ALLOCATE (ISG(MSA),ISPF(MSA),IWTH(MSA),JCN(MSA),JCN0(MSA),JCN1(MSA),&
      JD(MSA),KC(MSA),KF1(MSA),KP1(MSA),KT(MSA),KTF(MSA),KTMX(MSA),KTT&
      (MSA),LM(MSA),LRD(MSA),LUN(MSA),LUNS(MSA),MXSR(MSA),NBCF(MSA),NBCT&
      (MSA),NBFF(MSA),NBFT(MSA),NBSA(MSA),NBSL(MSA),NBW(MSA),NDFA(MSA),&
      NII(MSA),NMW(MSA),NRO(MSA),NVCN(MSA),NWDA(MSA))
      ALLOCATE (IPSF(MPO),KPSN(MPO),JPC(MPS),KPC(MPS),NPC(MPS),IHX(MHX)) 
	  ALLOCATE (IHU(MNC,MSA),IYH(MNC,MSA),JE(MNC,MSA),JP(MNC,MSA),&                       
      JPL(MNC,MSA),KGO(MNC,MSA),NCR(MNC,MSA),NHU(MNC,MSA),NYLN(MNC,MSA))             
	  ALLOCATE (NCP(MRO,MSA),NFRT(MRO,MSA),NPST(MRO,MSA),NTL(MRO,MSA))                    
	  ALLOCATE (IDFA(MHD,MOW),IDFD(MHD,MOW),IDMU(MHD,MOW),IGZO(MHD,MOW),&                 
      IGZX(MHD,MOW),IHBS(MHD,MOW),IYHO(MHD,MOW),LGIR(MHD,MOW),NBHS(MHD,&             
      MOW),NCOW(MHD,MOW),NGZA(MHD,MOW),NHBS(MHD,MOW),NYHO(MHD,MOW))                  
	  ALLOCATE (IDSL(MSA,MSA),IDSS(MSA,MSA),IDOW(MSA,MOW),&                 
      IHT(MNT,MSA),KOMP(MNT,MSA),IFED(MHD,MSA),NGZ(MHD,MSA),LID(ML1,MSA)&            
      ,LORG(MSL,MSA),IFLO(MSL,MSA),IHRL(12,MSA),IDFT(5,MSA),IDF0(5,MSA))                          
	  ALLOCATE (ITL(MRO,MNT,MSA),JH(MRO,MNT,MSA),KDT(12,MNC,MSA),LFT(MRO&                 
      ,MNT,MSA),LT(MRO,MNT,MSA),LYR(MRO,MNT,MSA),LPC(MRO,MPS,MSA),LY(MRO&            
      ,MNC,MSA),IHDT(MBS,MHD,MOW),NGIX(MSA,MHD,MOW),NBSX(MBS,MHD,MOW))                                 
      ALLOCATE (OSAA(MOW),OWSA(MOW),PKRZ(MSL),RNMN(MSL),RSPC(MSL),SSO3&
      (MSL),UN(MSL),UP(MSL),UW(MSL))
      ALLOCATE (FCST(MFT),FN(MFT),FNMA(MFT),FNMN(MFT),FNO(MFT),FOC(MFT),&
      FP(MFT),FPO(MFT),FSLT(MFT))                                                 
	  ALLOCATE (PCST(MPS),PHLF(MPS),PHLS(MPS),PKOC(MPS),PLCH(MPS),&                       
      PSOL(MPS),PWOF(MPS),SSPS(MPS))                                                 
      ALLOCATE (COOP(MNT),COTL(MNT),DKH(MNT),DKI(MNT),EFM(MNT),EMX&                  
      (MNT),FPOP(MNT),FRCP(MNT),FULU(MNT),HE(MNT),HMO(MNT),ORHI(MNT),RHT&            
      (MNT),RIN(MNT),RR(MNT),STIR(MNT),TIL(MNT),TLD(MNT))                                      
	  ALLOCATE (AEP(MNC),ALT(MNC),CAF(MNC),CNY(MNC),CSTS(MNC),CPY(MNC),&                  
      DDM(MNC),DLAI(MNC),DMLA(MNC),DMLX(MNC),EP(MNC),FLT(MNC),FTO(MNC),&
      GMHU(MNC),GSI(MNC),HI(MNC),HMX(MNC),PHUX(MNC),PLAX(MNC),POPX(MNC),&
      PRYF(MNC),PRYG(MNC),PST(MNC),RBMD(MNC),RDMX(MNC),RLAD(MNC),SDW(MNC),&
      TBSC(MNC),TCPA(MNC),TCPY(MNC),TOPC(MNC),VPD2(MNC),VPTH(MNC),WA(MNC),&
      WAVP(MNC),WCY(MNC),WSYF(MNC),WXYF(MNC),XDLAI(MNC),XMTU(MNC),YLD(MNC))                    
	  ALLOCATE (ABD(MSA),AFLG(MSA),AGPM(MSA),ALGI(MSA),ALQ(MSA),ARMN&                     
      (MSA),ARMX(MSA),ARSD(MSA),BA1(MSA),BA2(MSA),BCOF(MSA),BCV(MSA),&
      BFFL(MSA),BFSN(MSA),BFT(MSA),BIG(MSA),BIR(MSA),BR1(MSA),BR2(MSA),&
      BTC(MSA),BTCX(MSA),BTN(MSA),BTNX(MSA),BTP(MSA),BTPX(MSA),BV1(MSA),&
      BV2(MSA),BVIR(MSA),CFNP(MSA),CHL(MSA),CHN(MSA),CHS(MSA),CHXA(MSA),&
      CHXP(MSA),CLG(MSA),CN0(MSA),CN2(MSA),CNSX(MSA))
      ALLOCATE (COST(MSA),COWW(MSA),CST1(MSA),CV(MSA),CVF(MSA),CVP(MSA),&
      CVRS(MSA),CYAV(MSA),CYMX(MSA),CYSD(MSA),DALG(MSA),DDLG(MSA),DEPC&
      (MSA),DHT(MSA),DKIN(MSA),DKHL(MSA),DRT(MSA),DST0(MSA),DWOC(MSA),&
      EFI(MSA),EK(MSA),EM10(MSA),FBM(MSA),FCMN(MSA),FCMP(MSA),FDSF(MSA),&
      FFC(MSA),FFPQ(MSA),FGC(MSA),FGSL(MSA),FHP(MSA),FIRG(MSA),FPF(MSA),&
      FSFN(MSA),FSFP(MSA),GMA(MSA),GRDL(MSA),GWSN(MSA),GWST(MSA),GWMX(MSA),&
      HCLD(MSA),HCLN(MSA),HLMN(MSA),HR0(MSA),HSM(MSA),OCPD(MSA))
      ALLOCATE (OMAP(MSA),ORSD(MSA),PAW(MSA),PCOF(MSA),PDAW(MSA),PDPL&
      (MSA),PDPL0(MSA),PDPLC(MSA),PDSW(MSA),PEC(MSA),PMX(MSA),PM10(MSA),&
      PRSD(MSA),PSTF(MSA),QCAP(MSA),QRBQ(MSA),QRQB(MSA),RCBW(MSA),RCF&
      (MSA),RCHC(MSA),RCHD(MSA),RCHK(MSA),RCHL(MSA),RCHN(MSA),RCHS(MSA),&
      RCHX(MSA),RCSS(MSA),RCTW(MSA),REPI(MSA),RFPK(MSA),RFPL(MSA),RFPS&
      (MSA),RFPW(MSA),RFPX(MSA),RFTT(MSA),RFV(MSA),RFV0(MSA),RHD(MSA),&
      RHTT(MSA))
      ALLOCATE (RINT(MSA),RLF(MSA),RMXS(MSA),ROSP(MSA),RRUF(MSA),RSAE&
      (MSA),RSAP(MSA),RSBD(MSA),RSDP(MSA),RSEE(MSA),RSEP(MSA),RSF(MSA),&
      RSHC(MSA),RSK(MSA),RSLK(MSA),RSOC(MSA),RSON(MSA),RSOP(MSA),RSO3&
      (MSA),RSRR(MSA),RSSA(MSA),RSSP(MSA),RSV(MSA),RSVB(MSA),RSVE(MSA),&
      RSVF(MSA),RSVP(MSA),RSYB(MSA),RSYF(MSA),RSYN(MSA),RSYS(MSA),RVE0&
      (MSA),RVP0(MSA),RZ(MSA),RZSW(MSA),SALB(MSA),SAMA(MSA),SATK(MSA),&
      SCI(MSA),SDVR(MSA),SLF(MSA),SLT0(MSA),SMAS(MSA),SMEO(MSA),SMFN(MSA),&
      SMFU(MSA),SMLA(MSA),SMMU(MSA),SMNU(MSA),SMPL(MSA),SMPQ(MSA),SMPY&
      (MSA),SMRF(MSA),SMSS(MSA),SMST(MSA),SMWS(MSA),SMX(MSA),SMY1(MSA),&
      SMY2(MSA),SNO(MSA),SOLQ(MSA),SPLG(MSA),SRAD(MSA),SRSD(MSA),&
      SSFI(MSA),SSIN(MSA),SSW(MSA))                                                  
      ALLOCATE (ST0(MSA),STDO(MSA),STDON(MSA),STDOP(MSA),STKR(MSA),STLT&
      (MSA),STP(MSA),SW(MSA),SWB(MSA),SWLT(MSA),S3(MSA),TAGP(MSA),TCC(MSA),&               
      TCS(MSA),TFLG(MSA),THK(MSA),TILG(MSA),TKR(MSA),TLMF(MSA),TMN(MSA),&            
      TMX(MSA),TNOR(MSA),TOC(MSA),TRSD(MSA),TSLA(MSA),TSMY(MSA),TSNO(MSA),&
      TVGF(MSA),TYN(MSA),TYP(MSA),U10(MSA),UB1(MSA),UOB(MSA),UPSX(MSA),&
      URBF(MSA),USL(MSA),VAC(MSA),VALF1(MSA),VAP(MSA),VCHA(MSA),VCHB(MSA),&
      VFPA(MSA),VFPB(MSA),VIMX(MSA),VIRT(MSA),VIR0(MSA),VLG(MSA),VLGB&
      (MSA),VLGI(MSA),VLGM(MSA),VLGN(MSA),VPU(MSA),VRSE(MSA),VSLT(MSA),&
      WDRM(MSA))                                                  
	  ALLOCATE (WK(MSA),WS(MSA),WSA(MSA),WSX(MSA),WTMB(MSA),WTBL(MSA),&                   
      WTMN(MSA),WTMU(MSA),WTMX(MSA),XCT(MSA),XHSM(MSA),XIDK(MSA),XIDS&               
      (MSA),XMAP(MSA),XNS(MSA),XRFI(MSA),XTP1(MSA),YCT(MSA),YLC(MSA),&
      YLS(MSA),YTN(MSA),ZBMC(MSA),ZBMN(MSA),ZCO(MSA),ZCOB(MSA),ZFOP(MSA),&
      ZHPC(MSA),ZHPN(MSA),ZHSC(MSA),ZHSN(MSA),ZLM(MSA),ZLMC(MSA),ZLMN(MSA),&
      ZLS(MSA),ZLSC(MSA),ZLSL(MSA),ZLSLC(MSA),ZLSLNC(MSA),ZLSN(MSA),ZNMA&
      (MSA),ZNMN(MSA),ZNMU(MSA),ZNOA(MSA),ZNOS(MSA),ZNOU(MSA),ZOC(MSA),&
      ZON(MSA),ZPMA(MSA),ZPML(MSA),ZPMS(MSA),ZPMU(MSA),ZPO(MSA),ZPOU(MSA),&
      ZSLT(MSA),ZTP(MSA))                                             
	  ALLOCATE (CPVH(MHY),DPMT(MHY),DRAV(MHY),ERAV(MHY),HYDV(MHY),RCTC&                   
      (MHY),PRAV(MHY),PRB(MHY),PSZM(MHY),QC(MHY),QDR(MHY),QDRN(MHY),&             
      QN(MHY),QP(MHY),QPR(MHY),QPU(MHY),QRF(MHY),QRFN(MHY),QRP(MHY),&
      QURB(MHY),QVOL(MHY),RQRB(MHY),RSFN(MHY),RSSF(MHY),RWSA(MHY),&
      SHYD(MHY),SQVL(MHY),SST(MHY),STY(MHY),TC(MHY),TCAV(MHY),TCMN(MHY),&
      TCMX(MHY),TSFN(MHY),WYLD(MHY),YC(MHY),YCOU(MHY),YCWN(MHY),YMNU&
      (MHY),YN(MHY),YNOU(MHY),YNWN(MHY),YP(MHY),YPOU(MHY),YPWN(MHY),&
      YW(MHY))
      ALLOCATE (PSO3(MPO),PSON(MPO),PSOP(MPO),PSOQ(MPO),PSOY(MPO),PQPS&
      (MPO),PSSP(MPO),PYPS(MPO))
      ALLOCATE (QGA(MHP),RFDT(MHP))
      ALLOCATE (EO5(30,MSA),RF5(30,MSA),SCFS(30,MSA),XMS(30,MSA),ASW(12,&
      MSA),CX(12,MSA),QIN(12,MSA),SET(12,MSA),SRD(12,MSA),SRMX(12,MSA),&
      TAMX(12,MSA),TCN(12,MSA),TCVF(12,MSA),TEI(12,MSA),TET(12,MSA),THRL&
      (12,MSA),TQ(12,MSA),TR(12,MSA),TRHT(12,MSA),TSN(12,MSA),TSR(12,MSA),&
      TSY(12,MSA),TXMX(12,MSA),TXMN(12,MSA),TQN(12,MSA),TQP(12,MSA),TQPU&
      (12,MSA),TYON(12,MSA),TYTP(12,MSA),TYW(12,MSA),CNSC(2,MSA))
      ALLOCATE (ALS(MSL,MSA),BD(MSL,MSA),BDD(MSL,MSA),BDM(MSL,MSA),BDP&              
      (MSL,MSA),BK(MSL,MSA),CAC(MSL,MSA),CBN(MSL,MSA),CEC(MSL,MSA),CLA&              
      (MSL,MSA),CNDS(MSL,MSA),CNRT(MSL,MSA),CPRH(MSL,MSA),CPRV(MSL,MSA),&
      DHN(MSL,MSA),ECND(MSL,MSA),EXCK(MSL,MSA),FC(MSL,MSA),FOP(MSL,MSA),&
      HCL(MSL,MSA),PH(MSL,MSA),PO(MSL,MSA),PSP(MSL,MSA),ROK(MSL,MSA),RSD&
      (MSL,MSA),RSDM(MSL,MSA),SAN(MSL,MSA),SATC(MSL,MSA),SEV(MSL,MSA),&
      SIL(MSL,MSA),SMB(MSL,MSA))            
      ALLOCATE (ST(MSL,MSA),STFR(MSL,MSA),STMP(MSL,MSA),S15(MSL,MSA),WBMC&
      (MSL,MSA),WBMN(MSL,MSA),WCMU(MSL,MSA),WCOU(MSL,MSA),WHPC(MSL,MSA),&
      WHPN(MSL,MSA),WHSC(MSL,MSA),WHSN(MSL,MSA),WNMA(MSL,MSA),WLM(MSL,&
      MSA),WLMC(MSL,MSA),WLMN(MSL,MSA),WLS(MSL,MSA),WLSC(MSL,MSA),WLSL&
      (MSL,MSA),WLSLC(MSL,MSA),WLSLNC(MSL,MSA),WLSN(MSL,MSA),WNMN(MSL,MSA)&
      ,WNMU(MSL,MSA),WNOU(MSL,MSA),WOC(MSL,MSA),WON(MSL,MSA),WPMA(MSL,MSA),&              
      WPML(MSL,MSA),WPMS(MSL,MSA),WPMU(MSL,MSA),WPO(MSL,MSA),WPOU(MSL,&              
      MSA),WSLT(MSL,MSA),WT(MSL,MSA),Z(MSL,MSA))                                                   
      ALLOCATE (ACET(MNC,MSA),AJHI(MNC,MSA),AWC(MNC,MSA),CAW(MNC,MSA),&              
      CPHT(MNC,MSA),CSTF(MNC,MSA),DM(MNC,MSA),DMF(MNC,MSA),DM1(MNC,MSA),&            
      ETG(MNC,MSA),FRTN(MNC,MSA),FRTP(MNC,MSA),HU(MNC,MSA),HUF(MNC,MSA),&
      HUI(MNC,MSA),PPL0(MNC,MSA),RD(MNC,MSA),RDF(MNC,MSA),REG(MNC,MSA),&
      RW(MNC,MSA),SLAI(MNC,MSA),SLA0(MNC,MSA),SRA(MNC,MSA),STD(MNC,MSA),&
      STDL(MNC,MSA),STDN(MNC,MSA))               
      ALLOCATE (STDP(MNC,MSA),STL(MNC,MSA),SWH(MNC,MSA),SWP(MNC,MSA),TCAW&
      (MNC,MSA),TDM(MNC,MSA),TETG(MNC,MSA),TFTN(MNC,MSA),TFTP(MNC,MSA),&
      THU(MNC,MSA),TRA(MNC,MSA),TRD(MNC,MSA),TYL1(MNC,MSA),TYL2(MNC,MSA),&
      TYLN(MNC,MSA),TYLP(MNC,MSA),UNA(MNC,MSA),UN1(MNC,MSA),UP1(MNC,MSA),&
      VIR(MNC,MSA),WLV(MNC,MSA),XLAI(MNC,MSA),XDLA0(MNC,MSA),YLD1(MNC,&
      MSA),YLD2(MNC,MSA),YLNF(MNC,MSA),YLPF(MNC,MSA))     
      ALLOCATE (YHY(MHP,MHY),SMH(NSH,MHY),SMYH(NSH,MHY),&               
      VARH(NSH,MHY),QPST(MPS,MHY),RSPS(MPS,MHY),TSPS(MPS,MHY),YPST(MPS,&
      MHY),SRCH(26,MHY),CPFH(MSL,MHY),QSF(MSL,MHY),SSF(MSL,MHY),SM(NSM,&
      MSA),SMY(NSM,MSA),VAR(NSM,MSA),CTSA(100,MSA),VQ(90,MHY),VY(90,MHY),&
      GWPS(MPS,MSA),PFOL(MPS,MSA),ANA(MRO,MSA),FNMX(MRO,MSA),GCOW(MHD,MSA),&
      GZLM(MHD,MSA),DUMP(MHD,MOW),FFED(MHD,MOW),GZRT(MHD,MOW),VURN(MHD,MOW),&
      PPX(13,MSA),YSD(8,MHY),PCT(5,MHY),PCTH(5,MHY),SQB(5,MHY),SYB(5,MHY),&
      FNP(5,MSA),SMYRP(5,MPS),BN(4,MNC),BP(4,MNC),BLG(3,MNC),BWN(3,MNC),&
      DLAP(2,MNC),FRST(2,MNC),PPCF(2,MNC),PPLP(2,MNC),RWPC(2,MNC),STX&
      (2,MNC),WAC2(2,MNC))    
 	  ALLOCATE (SMMH(NSH,12,MHY),PVQ(MPS,90,MHY),PVY(MPS,90,MHY),PSTE&                   
      (MRO,MPS,MSA),PSTR(MRO,MPS,MSA),PSSF(MPS,MSL,MSA),PSTZ(MPS,MSL,&               
      MSA),CND(MRO,MNT,MSA),SMM(NSM,12,MSA),STV(17,12,MSA),SMS(11,ML1,&               
      MSA),SOL(7,MSL,MSA),FIRX(MRO,MNT,MSA),QIR(MRO,MNT,MSA),RSTK(MRO,&
      MNT,MSA),TIR(MRO,MNT,MSA),VIRR(MRO,MNT,MSA),FDP(MRO,MNT,MSA),WFA&
      (MRO,MFT,MSA),PHU(MNC,MRO,MSA),POP(MNC,MRO,MSA),PPLA(MNC,MRO,MSA),&
      VARC(16,MNC,MSA),SMAP(13,MPS,MHY),SMYP(13,MPS,MHY),VARP(12,MPS,MHY),&
      SMMRP(5,MPS,12),SMRP(5,MPS,12),SPQ(5,20,MHY),SPQC(5,20,MHY),SPY(5,&
      20,MHY),TSFC(6,MNC,MSA),SOIL(17,MSL,MSA),HUSC(MRO,MNT,MSA),RWT(MSL,&
      MNC,MSA),STDA(4,MNC,MSA),SFCP(6,MNC,MSA),SFMO(5,MNC,MSA),&
      XZP(13,ML1,MSA),QHY(NPD,MHY,MHX))
      ALLOCATE (SMMP(20,MPS,13,MHY),SMMC(15,MNC,12,MSA))
!     ALLOCATE (APQ(5,20,100,MHY),APQC(5,20,100,MHY),APY(5,20,100,MHY)&              
!    &,AQB(5,20,100,MHY),AYB(5,20,100,MHY),SMMP(20,MPS,12,MSA),SMMC(15,&             
!    &MNC,12,MSA))                                                                                            
      RETURN                                                                         
      END                                                                            
                                                                                     
