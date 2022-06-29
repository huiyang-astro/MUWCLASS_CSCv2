# MUWCLASS_CSCv2
 
## Classifying Unidentified X-ray Sources in the Chandra Source Catalog Using a Multi-wavelength Machine Learning Approach
### Hui Yang1, Jeremy Hare2, Oleg Kargaltsev1, Igor Volkov1, Steven Chen1, Blagoy Rangelov3
### 1 The George Washington University 2 NASA GSFC 3 Texas State University

## CHECK our paper at https://arxiv.org/abs/2206.13656

### contact huiyang@gwu.edu if you have any questions

This github repo provides the MUltiWavelength Machine Learning CLASSification Pipeline (MUWCLASS) and the classification results on the Chandra Source Catalog v2 (CSCv2).

The main components of this github repo are

codes/demo
- This is a notebook of a demonstration of classifying CSCv2 sources using MUWCLASS with CSCv2 and multiwavelength data

codes/builtTD
- This is a notebook that rebuild the training dataset used for the MUWCLASS

files/HESS_results
- This includes the classification results of X-ray sources within the extent of unidentified HESS sources with detection significance S/N>6 (files/HESS_results/signif_6) and S/N>3 (files/HESS_results/signif_3)

-- 

* This notebook was run in CIAO 4.14 with Python 3.9 
* run the follow code to create a new conda environment ciao-4.14-muwclass; if you already have ciao-4.14 installed with Python 3.9, you can use your own conda environment with additional Python packages installed from below
* conda create -n ciao-4.14-muwclass -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge ciao sherpa ds9 ciao-contrib caldb_main marx python=3.9

* run 'bash install-packages.sh' under ciao-4.14-muwclass environment to install all required packages 

* then, make sure to enable widgetsnbextension and ipyaladin, run 
* jupyter nbextension enable --py widgetsnbextension
* jupyter nbextension enable --py --sys-prefix ipyaladin
- on your terminal 

* You might also need to manually register the existing ds9 with the xapns name server by selecting the ds9 File->XPA->Connect menu option so your ds9 will be fully accessible to pyds9. 

