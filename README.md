# SPACE-MS overview 

SPACE-MS is a repository of command line tools that enable spatial analysis of tissue damage distribution in multiple sclerosis.

# SPACE-MS dependencies
To run SPACE-MS you need Python 3 (for example, via [Anaconda](http://www.anaconda.com/distribution)). You will also need:
* [NumPy](http://numpy.org)
* [Nibabel](http://nipy.org/nibabel)
* [PyKrige](http://pypi.org/project/PyKrige)

# SPACE-MS download
To get SPACE-MS (Linux or MacOs):

1. Open a terminal and go to the folder where you want to store SPACE-MS;
2. Clone the repository:
```
git clone https://github.com/carmentur/SPACE-MS.git 
```
3. You will now have the SPACE-MS tools ready for you in the `SPACE-MS/spacetools` folder. 

4. To check how to use a tool, simply print its manual by typing in your terminal:
```
python </PATH/TO/TOOL> --help
```

# SPACE-MS: list of tools

The following tools are available:
* `run_3dvario.py`: to perform 3D variogram analysis of a 3D NIFTI file storing a lesion mask;
* ... more tools coming soon - stay tuned!


# Do you use SPACE-MS?
If you use SPACE-MS tools, please cite us:

"Linking macrostructural and microstructural damage in early MS: a geostatistical and diffusion MRI study". Carmen Tur, Robert Marschallinger, Ferran Prados, Sara Collorone, Daniel R Altmann, Sébastien Ourselin, Claudia A. M. Gandini Wheeler-Kingshott and Olga Ciccarelli; Proceedings of the 2018 meeting of the ISMRM, p. 0090.

# License
SPACE-MS tools are distributed under the BSD 2-Clause License, Copyright (c) 2020, University College London. All rights reserved (license at this [link](https://github.com/carmentur/SPACE-MS/blob/master/LICENSE.txt)).

# Acknowledgements
2015 European Committee for Treatment and Research in Multiple Sclerosis (ECTRIMS) post-doctoral research fellowship; UK Multiple Sclerosis Society (grants 892/08 and 77/2017). The help of Dr [Francesco Grussu](http://github.com/fragrussu) is also acknowledged. 
