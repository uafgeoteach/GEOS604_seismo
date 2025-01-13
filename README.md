# GEOS 604 &mdash; Seismology   

**Course:** GEOS604   
**University:** University of Alaska Fairbanks    
**Instructor:** Bryant Chow  
**Website:** https://bryantchow.com/teaching/geos604

## Setup

To do the homework assignments in this repository, you will need a [Conda](https://anaconda.org/anaconda/conda) environment. I recommend installing [Miniconda](https://docs.anaconda.com/miniconda/install/#quick-command-line-install)
Once you have installed Miniconda, run the following commands to clone the class repository and install the class environment.
```bash
cd <PATH_TO_WORKING_DIRECTORY>
git clone git@github.com:uafgeoteach/GEOS604_seismo.git
cd GEOS604_seismo
conda env create -f environment.yml
conda activate geos604
```

To run the notebooks, we will be using [Jupyter Lab](https://jupyter.org/) which has been installed in the Conda environment. To start Jupyter lab you simply need to type
```bash
jupyter lab
```

All homework assignments are located in `homeworks/` as either Markdown files or Jupyter notebooks. Markdown files are automatically rendered when viewed in GitHub. 
