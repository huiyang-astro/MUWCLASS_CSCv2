# files

This directory contains 

HESS_results
- This includes the classification results of X-ray sources within the extent of unidentified HESS sources with detection significance S/N>6 (files/HESS_results/signif_6) and S/N>3 (files/HESS_results/signif_3)

MRT
- The three machine readable tables for the properties and the classification results of GCS, TD and HESS field sources. 

CSC_TD_MW_remove.csv
- The current training dataset used for the MUWCLASS

You also need to run the following python codes to clean the TD before using it.
* TD = pd.read_csv('CSC_TD_MW_remove.csv')
* TD['Class'] = TD['Class'].replace({'NS_BIN':'LMXB'})
* TD = prepare_cols(TD, cp_thres=0, TD=True, NS_MWdrop=False, STAR_classremove=['HM-STAR','LM-STAR','YSO']) # some filtering 

LOO_classes.csv
- The classification results of the TD Using LOOCV

feature_importance.csv
- The feature importance of the all possible features

figs/
- some figures

