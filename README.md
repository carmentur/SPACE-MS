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
```
./MyRelax/myrelax
```
4. To check how to use a tool, simply print its manual by typing in your terminal:
```
python </PATH/TO/TOOL> --help
```

# SPACE-MS: list of tools

The following tools are available:
* `run_3dvario.py`: to perform 3D variogram analysis of a 3D NIFTI file storing a lesion mask;
* ... more tools coming soon - stay tuned!

