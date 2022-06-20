# Build_TD

This notebook is used to reconstruct the current training dataset (files/CSC_TD_MW_remove.csv) used for the MUWCLASS. The reconstructed TD will be saved in codes/buildTD/data/CSC_TD_MW_remove.csv. It is noted that you may get a slightly different TD using a different computer or different versions of python or python packages. 

* Please follow the installation instructions on the main page of this github repo to install all required packages. 

You also need to run the following python codes to clean the TD before using it.
* TD = pd.read_csv('CSC_TD_MW_remove.csv')
* TD['Class'] = TD['Class'].replace({'NS_BIN':'LMXB'})
* TD = prepare_cols(TD, cp_thres=0, TD=True, NS_MWdrop=False, STAR_classremove=['HM-STAR','LM-STAR','YSO']) # some filtering 

