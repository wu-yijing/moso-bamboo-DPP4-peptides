@echo off
REM 批量 Vina 对接 - 用 Docker 或本地 Vina
set REC=1WCY_receptor.pdbqt
set BOX=moso_box.txt
echo --- Docking CPR ---
vina --receptor %REC% --ligand moso_ligands\00_CPR.pdbqt --config %BOX% --out moso_ligands\dock_00_CPR.pdbqt --cpu 4
echo --- Docking CPV ---
vina --receptor %REC% --ligand moso_ligands\01_CPV.pdbqt --config %BOX% --out moso_ligands\dock_01_CPV.pdbqt --cpu 4
echo --- Docking IAP ---
vina --receptor %REC% --ligand moso_ligands\02_IAP.pdbqt --config %BOX% --out moso_ligands\dock_02_IAP.pdbqt --cpu 4
echo --- Docking IPA ---
vina --receptor %REC% --ligand moso_ligands\03_IPA.pdbqt --config %BOX% --out moso_ligands\dock_03_IPA.pdbqt --cpu 4
echo --- Docking IPE ---
vina --receptor %REC% --ligand moso_ligands\04_IPE.pdbqt --config %BOX% --out moso_ligands\dock_04_IPE.pdbqt --cpu 4
echo --- Docking IPI ---
vina --receptor %REC% --ligand moso_ligands\05_IPI.pdbqt --config %BOX% --out moso_ligands\dock_05_IPI.pdbqt --cpu 4
echo --- Docking IPK ---
vina --receptor %REC% --ligand moso_ligands\06_IPK.pdbqt --config %BOX% --out moso_ligands\dock_06_IPK.pdbqt --cpu 4
echo --- Docking IPR ---
vina --receptor %REC% --ligand moso_ligands\07_IPR.pdbqt --config %BOX% --out moso_ligands\dock_07_IPR.pdbqt --cpu 4
echo --- Docking IPS ---
vina --receptor %REC% --ligand moso_ligands\08_IPS.pdbqt --config %BOX% --out moso_ligands\dock_08_IPS.pdbqt --cpu 4
echo --- Docking IPV ---
vina --receptor %REC% --ligand moso_ligands\09_IPV.pdbqt --config %BOX% --out moso_ligands\dock_09_IPV.pdbqt --cpu 4
echo --- Docking LAP ---
vina --receptor %REC% --ligand moso_ligands\10_LAP.pdbqt --config %BOX% --out moso_ligands\dock_10_LAP.pdbqt --cpu 4
echo --- Docking LPA ---
vina --receptor %REC% --ligand moso_ligands\11_LPA.pdbqt --config %BOX% --out moso_ligands\dock_11_LPA.pdbqt --cpu 4
echo --- Docking LPC ---
vina --receptor %REC% --ligand moso_ligands\12_LPC.pdbqt --config %BOX% --out moso_ligands\dock_12_LPC.pdbqt --cpu 4
echo --- Docking LPD ---
vina --receptor %REC% --ligand moso_ligands\13_LPD.pdbqt --config %BOX% --out moso_ligands\dock_13_LPD.pdbqt --cpu 4
echo --- Docking LPE ---
vina --receptor %REC% --ligand moso_ligands\14_LPE.pdbqt --config %BOX% --out moso_ligands\dock_14_LPE.pdbqt --cpu 4
echo --- Docking LPG ---
vina --receptor %REC% --ligand moso_ligands\15_LPG.pdbqt --config %BOX% --out moso_ligands\dock_15_LPG.pdbqt --cpu 4
echo --- Docking LPH ---
vina --receptor %REC% --ligand moso_ligands\16_LPH.pdbqt --config %BOX% --out moso_ligands\dock_16_LPH.pdbqt --cpu 4
echo --- Docking LPI ---
vina --receptor %REC% --ligand moso_ligands\17_LPI.pdbqt --config %BOX% --out moso_ligands\dock_17_LPI.pdbqt --cpu 4
echo --- Docking LPK ---
vina --receptor %REC% --ligand moso_ligands\18_LPK.pdbqt --config %BOX% --out moso_ligands\dock_18_LPK.pdbqt --cpu 4
echo --- Docking LPN ---
vina --receptor %REC% --ligand moso_ligands\19_LPN.pdbqt --config %BOX% --out moso_ligands\dock_19_LPN.pdbqt --cpu 4
echo --- Docking LPP ---
vina --receptor %REC% --ligand moso_ligands\20_LPP.pdbqt --config %BOX% --out moso_ligands\dock_20_LPP.pdbqt --cpu 4
echo --- Docking LPQ ---
vina --receptor %REC% --ligand moso_ligands\21_LPQ.pdbqt --config %BOX% --out moso_ligands\dock_21_LPQ.pdbqt --cpu 4
echo --- Docking LPR ---
vina --receptor %REC% --ligand moso_ligands\22_LPR.pdbqt --config %BOX% --out moso_ligands\dock_22_LPR.pdbqt --cpu 4
echo --- Docking LPS ---
vina --receptor %REC% --ligand moso_ligands\23_LPS.pdbqt --config %BOX% --out moso_ligands\dock_23_LPS.pdbqt --cpu 4
echo --- Docking LPT ---
vina --receptor %REC% --ligand moso_ligands\24_LPT.pdbqt --config %BOX% --out moso_ligands\dock_24_LPT.pdbqt --cpu 4
echo --- Docking LPV ---
vina --receptor %REC% --ligand moso_ligands\25_LPV.pdbqt --config %BOX% --out moso_ligands\dock_25_LPV.pdbqt --cpu 4
echo --- Docking VPA ---
vina --receptor %REC% --ligand moso_ligands\26_VPA.pdbqt --config %BOX% --out moso_ligands\dock_26_VPA.pdbqt --cpu 4
echo --- Docking VPG ---
vina --receptor %REC% --ligand moso_ligands\27_VPG.pdbqt --config %BOX% --out moso_ligands\dock_27_VPG.pdbqt --cpu 4
echo --- Docking VPI ---
vina --receptor %REC% --ligand moso_ligands\28_VPI.pdbqt --config %BOX% --out moso_ligands\dock_28_VPI.pdbqt --cpu 4
echo --- Docking VPK ---
vina --receptor %REC% --ligand moso_ligands\29_VPK.pdbqt --config %BOX% --out moso_ligands\dock_29_VPK.pdbqt --cpu 4
echo --- Docking VPR ---
vina --receptor %REC% --ligand moso_ligands\30_VPR.pdbqt --config %BOX% --out moso_ligands\dock_30_VPR.pdbqt --cpu 4
echo --- Docking VPT ---
vina --receptor %REC% --ligand moso_ligands\31_VPT.pdbqt --config %BOX% --out moso_ligands\dock_31_VPT.pdbqt --cpu 4
echo --- Docking VPV ---
vina --receptor %REC% --ligand moso_ligands\32_VPV.pdbqt --config %BOX% --out moso_ligands\dock_32_VPV.pdbqt --cpu 4
echo --- Docking APPR ---
vina --receptor %REC% --ligand moso_ligands\33_APPR.pdbqt --config %BOX% --out moso_ligands\dock_33_APPR.pdbqt --cpu 4
echo --- Docking CAPP ---
vina --receptor %REC% --ligand moso_ligands\34_CAPP.pdbqt --config %BOX% --out moso_ligands\dock_34_CAPP.pdbqt --cpu 4
echo --- Docking IPID ---
vina --receptor %REC% --ligand moso_ligands\35_IPID.pdbqt --config %BOX% --out moso_ligands\dock_35_IPID.pdbqt --cpu 4
echo --- Docking LPGP ---
vina --receptor %REC% --ligand moso_ligands\36_LPGP.pdbqt --config %BOX% --out moso_ligands\dock_36_LPGP.pdbqt --cpu 4
echo --- Docking LPPA ---
vina --receptor %REC% --ligand moso_ligands\37_LPPA.pdbqt --config %BOX% --out moso_ligands\dock_37_LPPA.pdbqt --cpu 4
echo --- Docking LPPG ---
vina --receptor %REC% --ligand moso_ligands\38_LPPG.pdbqt --config %BOX% --out moso_ligands\dock_38_LPPG.pdbqt --cpu 4
echo --- Docking LPPH ---
vina --receptor %REC% --ligand moso_ligands\39_LPPH.pdbqt --config %BOX% --out moso_ligands\dock_39_LPPH.pdbqt --cpu 4
echo --- Docking LPPQ ---
vina --receptor %REC% --ligand moso_ligands\40_LPPQ.pdbqt --config %BOX% --out moso_ligands\dock_40_LPPQ.pdbqt --cpu 4
echo --- Docking LPPV ---
vina --receptor %REC% --ligand moso_ligands\41_LPPV.pdbqt --config %BOX% --out moso_ligands\dock_41_LPPV.pdbqt --cpu 4
echo --- Docking LPSP ---
vina --receptor %REC% --ligand moso_ligands\42_LPSP.pdbqt --config %BOX% --out moso_ligands\dock_42_LPSP.pdbqt --cpu 4
echo --- Docking MPPA ---
vina --receptor %REC% --ligand moso_ligands\43_MPPA.pdbqt --config %BOX% --out moso_ligands\dock_43_MPPA.pdbqt --cpu 4
echo --- Docking VPPK ---
vina --receptor %REC% --ligand moso_ligands\44_VPPK.pdbqt --config %BOX% --out moso_ligands\dock_44_VPPK.pdbqt --cpu 4
echo --- Docking VPPN ---
vina --receptor %REC% --ligand moso_ligands\45_VPPN.pdbqt --config %BOX% --out moso_ligands\dock_45_VPPN.pdbqt --cpu 4
echo --- Docking APPSQ ---
vina --receptor %REC% --ligand moso_ligands\46_APPSQ.pdbqt --config %BOX% --out moso_ligands\dock_46_APPSQ.pdbqt --cpu 4
echo --- Docking APQIP ---
vina --receptor %REC% --ligand moso_ligands\47_APQIP.pdbqt --config %BOX% --out moso_ligands\dock_47_APQIP.pdbqt --cpu 4
echo --- Docking APSPE ---
vina --receptor %REC% --ligand moso_ligands\48_APSPE.pdbqt --config %BOX% --out moso_ligands\dock_48_APSPE.pdbqt --cpu 4
echo --- Docking APSQP ---
vina --receptor %REC% --ligand moso_ligands\49_APSQP.pdbqt --config %BOX% --out moso_ligands\dock_49_APSQP.pdbqt --cpu 4
echo --- Docking CPPSK ---
vina --receptor %REC% --ligand moso_ligands\50_CPPSK.pdbqt --config %BOX% --out moso_ligands\dock_50_CPPSK.pdbqt --cpu 4
echo --- Docking IPDAP ---
vina --receptor %REC% --ligand moso_ligands\51_IPDAP.pdbqt --config %BOX% --out moso_ligands\dock_51_IPDAP.pdbqt --cpu 4
echo --- Docking LAIPP ---
vina --receptor %REC% --ligand moso_ligands\52_LAIPP.pdbqt --config %BOX% --out moso_ligands\dock_52_LAIPP.pdbqt --cpu 4
echo --- Docking LAPSP ---
vina --receptor %REC% --ligand moso_ligands\53_LAPSP.pdbqt --config %BOX% --out moso_ligands\dock_53_LAPSP.pdbqt --cpu 4
echo --- Docking LPAQP ---
vina --receptor %REC% --ligand moso_ligands\54_LPAQP.pdbqt --config %BOX% --out moso_ligands\dock_54_LPAQP.pdbqt --cpu 4
echo --- Docking LPCPR ---
vina --receptor %REC% --ligand moso_ligands\55_LPCPR.pdbqt --config %BOX% --out moso_ligands\dock_55_LPCPR.pdbqt --cpu 4
echo --- Docking LPDDP ---
vina --receptor %REC% --ligand moso_ligands\56_LPDDP.pdbqt --config %BOX% --out moso_ligands\dock_56_LPDDP.pdbqt --cpu 4
echo --- Docking LPGDP ---
vina --receptor %REC% --ligand moso_ligands\57_LPGDP.pdbqt --config %BOX% --out moso_ligands\dock_57_LPGDP.pdbqt --cpu 4
echo --- Docking LPINP ---
vina --receptor %REC% --ligand moso_ligands\58_LPINP.pdbqt --config %BOX% --out moso_ligands\dock_58_LPINP.pdbqt --cpu 4
echo --- Docking LPPGP ---
vina --receptor %REC% --ligand moso_ligands\59_LPPGP.pdbqt --config %BOX% --out moso_ligands\dock_59_LPPGP.pdbqt --cpu 4
echo ALL DOCKING COMPLETE
pause
