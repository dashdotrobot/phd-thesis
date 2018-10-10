# PhD Dissertation: The Bicycle Wheel
__Author:__ Matthew Ford

__Institution:__ Northwestern University

__Department:__ Mechanical Engineering

_Copyright 2018 by Matthew Ford - All Rights Reserved_

## How to build the thesis document (thesis.pdf)
1. Install a LaTeX distribution such as [TeXLive](https://www.tug.org/texlive/)
2. Install [Inkscape](https://inkscape.org) and make sure that the `inkscape` command is available from the command line
3. Download or clone this repository to your local computer
4. Open your favorite command line in the same directory as `thesis.tex`
5. Run `latexmk -pdf -quiet --shell-escape thesis.tex`

Alternatively, use your favorite LaTeX distribution and build with bibtex support. Make sure that your build system allows the latex command `\write18`, which is required to convert the figures from .svg to .pdf_tex.

## How to reproduce the analysis
Most of the analysis in this thesis including some analytical derivations, plot generation, and numerical calculations is written in [Python](https://www.python.org/). Most of the code is in the form of [Jupyter notebooks](http://jupyter.org/) in the `code/` directory.

### 1. Install [Miniconda](https://conda.io/miniconda.html)

Make sure to download the Python 3.7 version.

### 2. Create a new environment with the necessary dependencies

Open an Anaconda Prompt or equivalent command line application in the directory containing `environment.yml` and execute the following command:

```
$ conda env create -f environment.yml
```

 This creates the `thesis` environment with all the necessary dependencies to run the code.

### 3. Activate the `thesis` environment and launch Jupyter Notebook

```
$ conda activate thesis
$ jupyter notebook
```

Now you can run the notebooks in `code/` in Jupyter!

### Note regarding ABAQUS simulations

Some of the Jupyter notebooks will generate ABAQUS input files under the `data/` directory. If you have ABAQUS installed, you can run the simulations by navigating to the directory containing the ABAQUS input files and executing all of the `_run_[i].bat` files (on Windows). After all the simulations have completed, run the `_postproc.bat` script to extract the relevant output. This will create a set of .csv files from each ABAQUS simulation. If you are running ABAQUS on Linux, you may need to convert the batch files to shell scripts.

If you are unable (or unwilling) to run the ABAQUS simulations, you can still run all of the plotting and analysis code in the `code/` directory. The .csv files produced by the postprocessing scripts are included in the repository.
