# MUWCLASS_CSCv2
 
## Classifying Unidentified X-ray Sources in the Chandra Source Catalog Using a Multi-wavelength Machine Learning Approach
### Hui Yang1, Jeremy Hare2, Oleg Kargaltsev1, Igor Volkov1, Steven Chen1, Blagoy Rangelov3
### 1 The George Washington University 2 NASA GSFC 3 Texas State University

## CHECK our paper at https://arxiv.org/abs/2206.13656

### contact huiyang@gwu.edu if you have any questions

This github repo provides the MUltiWavelength Machine Learning CLASSification Pipeline (MUWCLASS) and the classification results on the Chandra Source Catalog v2 (CSCv2).

The main components of this github repo are

files/{GCS.csv, GCS_MRT.dat, TD.csv, TD_MRT.dat, HESS.csv, HESS_MRT.dat}
- These include the CSV (comma-separated values) files of the good CSCv2 sample (GCS.csv), the training dataset (TD.csv) and the HESS field sources (HESS.csv) which are the full tables of Table 8, 9, and 10 in https://arxiv.org/abs/2206.13656, and their corresponding machine readable tables for the AAS Journals (GCS_MRT.dat, TD_MRT.dat, and HESS_MER.dat). 

codes/demo
- This is a notebook of a demonstration of classifying CSCv2 sources using MUWCLASS with CSCv2 and multiwavelength data

codes/builtTD
- This is a notebook that rebuild the training dataset used for the MUWCLASS

files/HESS_results
- This includes the classification results of X-ray sources within the extent of unidentified HESS sources with detection significance S/N>6 (files/HESS_results/signif_6) and S/N>3 (files/HESS_results/signif_3)

codes/Evaluation
- This is a notebook that evaluate the performance of the pipeline (calculating performance metrics, confusion matrices) for all the sources from the training dataset, or a subset of the sources that are missing optical/NIR/MIR counterparts. The users can adjust the filters to evaluate the performance any subset of the training dataset. 

files/{CSC_TD_MW_remove.csv, feature_importance.csv, tbabs.data}
- Some other CSV files including the raw training dataset with more properties (CSC_TD_MW_remove.csv), feature importance values and their uncertainties (feature_importance.csv), the photoelectric absorption cross-section file (tbabs.data).

-- Required Package Installation Instructions

* The notebooks were run in CIAO 4.14 with Python 3.9 

* If you are using macOS 11 (Big Sur) and macOS 12 (Monterey) on the new Apple M1 chip, you will need to use the Intel x86_64 edition of conda and Python for conda installation (see https://cxc.cfa.harvard.edu/ciao/download/platforms.html)

* XQuartz is required to install on macOS to run DS9.

* run the follow code to create a new conda environment ciao-4.14-muwclass; if you already have ciao-4.14 installed with Python 3.9, you can use your own conda environment with additional Python packages installed from below
* conda create -n ciao-4.14-muwclass -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge ciao sherpa ds9 ciao-contrib caldb_main marx python=3.9

* run 'bash install-packages.sh' under ciao-4.14-muwclass environment to install all required packages 

* then, make sure to enable widgetsnbextension and ipyaladin, run 
* jupyter nbextension enable --py widgetsnbextension
* jupyter nbextension enable --py --sys-prefix ipyaladin
- on your terminal 

* You might also need to manually register the existing ds9 with the xapns name server by selecting the ds9 File->XPA->Connect menu option so your ds9 will be fully accessible to pyds9. 

