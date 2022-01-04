# SPACE-MS overview 

SPACE-MS is a repository of python and matlab tools that enable spatial analysis of tissue damage distribution in multiple sclerosis.

# SPACE-MS dependencies
To run the python SPACE-MS tools you will need Python 3 (for example, via [Anaconda](http://www.anaconda.com/distribution)). You will also need:
* [NumPy](http://numpy.org)
* [Nibabel](http://nipy.org/nibabel)
* [Scikit-image](https://scikit-image.org)
* [PyKrige](http://pypi.org/project/PyKrige)

# SPACE-MS download
To get SPACE-MS (Linux or MacOs):

1. Open a terminal and go to the folder where you want to store SPACE-MS;
2. Clone the repository:
```
git clone https://github.com/carmentur/SPACE-MS.git 
```
3. You will now have the [SPACE-MS tools](https://github.com/carmentur/SPACE-MS/tree/master/spacetools) ready for you in the `SPACE-MS/spacetools` folder. 

4. To check how to use a tool, simply print its manual by typing in your terminal, e.g.
```
python run_sspace_ms.py --help
```
in python or 
```
help run_nci
```
in Matlab. 

5. You will also find a [data folder](https://github.com/carmentur/SPACE-MS/tree/master/data) storing some synthetic lesion masks on which you will be able to try out the SPACE-MS tools, alongside some tutorials on SPACE-MS (COMING SOON!). 



# SPACE-MS: list of tools

The following tools are available:
* [`run_3dvario.py`](https://github.com/carmentur/SPACE-MS/blob/master/spacetools/run_3dvario.py): to perform 3D variogram analysis of a 3D NIFTI file storing a lesion mask;
* [`run_sspace_ms.py`](https://github.com/carmentur/SPACE-MS/blob/master/spacetools/run_sspace_ms.py): to calculate SPACE-MS metrics on brain lesion masks as described in Carmen Tur et al, NeuroImage: Clinical 2022, 33: 102904, doi: [10.1016/j.nicl.2021.102904](https://doi.org/10.1016/j.nicl.2021.102904);
* [`run_nci.m`](https://github.com/carmentur/SPACE-MS/blob/master/spacetools/run_nci.m): to calculate the neruaxis caudality index (NCI) for brain lesions as described in Carmen Tur et al, NeuroImage: Clinical 2022, 33: 102904, doi: [10.1016/j.nicl.2021.102904](https://doi.org/10.1016/j.nicl.2021.102904) -- COMING SOON!


# Citations
If you use SPACE-MS tools, please cite our work as detailed below.

* If you use [`run_sspace_ms.py`](https://github.com/carmentur/SPACE-MS/blob/master/spacetools/run_sspace_ms.py) and/or [`run_nci.m`](https://github.com/carmentur/SPACE-MS/blob/master/spacetools/run_nci.m), please cite "Spatial patterns of brain lesions assessed through covariance estimations of lesional voxels in multiple sclerosis: the SPACE-MS technique". Carmen Tur et al, NeuroImage: Clinical 2022, 33: 102904, doi: [10.1016/j.nicl.2021.102904](https://doi.org/10.1016/j.nicl.2021.102904).

* If you use [`run_3dvario.py`](https://github.com/carmentur/SPACE-MS/blob/master/spacetools/run_3dvario.py), please cite "Linking macrostructural and microstructural damage in early MS: a geostatistical and diffusion MRI study". Carmen Tur, Robert Marschallinger, Ferran Prados, Sara Collorone, Daniel R Altmann, Sébastien Ourselin, Claudia A. M. Gandini Wheeler-Kingshott and Olga Ciccarelli; Proceedings of the 2018 meeting of the ISMRM, p. 0090.


# License
SPACE-MS tools are distributed under the BSD 2-Clause License, Copyright (c) 2020, 2021, 2022 University College London. All rights reserved (license at this [link](https://github.com/carmentur/SPACE-MS/blob/master/LICENSE.txt)).

**The use of SPACE-MS MUST also comply with the individual licenses of all of its dependencies.**

# Acknowledgements
This project received the support of a fellowship from ”la Caixa” Foundation (ID 100010434; the fellowship code is LCF/BQ/PI20/11760008) and from the 2021 Merck’s Award for the Investigation in Multiple Sclerosis (Spain). Funding from 2015 European Committee for Treatment and Research in Multiple Sclerosis (ECTRIMS) post-doctoral research fellowships and from the UK Multiple Sclerosis Society (grants 892/08 and 77/2017) are also acknowledged. The help of Dr [Francesco Grussu](http://fragrussu.github.io) is acknowledged. 
