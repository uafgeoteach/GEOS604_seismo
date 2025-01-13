# GEOS 604 &mdash; Seismology   

Course: GEOS604 
University: University of Alaska Fairbanks  
Instructor: Bryant Chow
Website: https://bryantchow.com/teaching/geos604

## Setup

To do the homework assignments in this repository, you will need a [Conda](https://anaconda.org/anaconda/conda) environment. I recommend installing [Miniconda](https://docs.anaconda.com/miniconda/install/#quick-command-line-install)
Once you have installed Miniconda, run the following commands to clone the class repository and install the class environment.
```bash
cd <PATH_TO_WORKING_DIRECTORY>
git clone git@github.com:uafgeoteach/GEOS604_seismo.git
cd GEOS604_seismo
conda env create -f environment.yml
```

All homework assignments are located in `homeworks/` as Markdown files. These should be rendered automatically when viewed in GitHub. I will also strive to convert Markdown files to PDFs and share those for those who may not have internet access, but please note that the markdown files will always be the most up-to-date source.

Students should fork this repository to their own GitHub account and complete and push Jupyter notebook homeworks in the `notebooks/` directory. See instructions inside HW1 Problem 0 to set up your working environment.
