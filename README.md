# A model for the role of STIM1 and STIM2 in Ca2+ oscillations in T-lymphocytes

## Description 


### The figure scripts

All the scripts used to generate the figures can be found in the PyScripts folder. PDF versions of the figures can be found in PyScripts/Figures.

To obtain a figure, simply run the corresponding script in python.

### The Python models

We use `Sympy` to integrate the ODE models. All the related functions used for time series integration and plotting can be found in the `STIMKO_Models.py` script. 

### Experimental data 

Snippet of the experimental data used in the figures can be found in individual csv files under `PyScripts/experimentalData`

### The AUTO scripts

Different versions of the model (WT, S1-KO, S2-KO) are containted in individual `.f90` files. Scripts with the  suffix `CDI` on the name also include the CDI regulatory mechanism. Scripts with the `DynamicCell` suffix on the name include both the CDI regulatory mechanism, the differential equation describing the changes in $K_1$ (STIM1's binding affinity to Ca2+), and the STIM-independent influx pathway ($\alpha_0$)

### The SymPy scripts

The code used to create the expressions for the activated forms $S_1$, $S_2$, $S_{21}$ and $S_{12}$ is found in the `PySym_STIMKOModel` and `PySym_STIMWTModel` scripts.  



## Requirements

The scripts are validated to run in python 3.8 using the following packages

- Matplotlib 3.4.0
- Numpy 1.19.1
- SciPy 1.9.0
- pandas 1.4.3

These version requirements have more to do with AUTO than the code itself (the interactive AUTO figure viewer opened by the `@pp` command does not seem to work with newer versions of matplotlib). A conda environment file (ca_env.yml) is provided in this repository with all the packages needed to run these scripts.

Make sure you have conda installed before importing the environment. The conda environment can be installed by typing:

```bash
conda env create -f ca_env.yml
```

One still needs to install AUTO to be able to run the continuation scripts, however.

## Installing AUTO

First, clone the AUTO repository. Open a terminal in the folder you wish to put the AUTO files in and type:

```bash
git clone https://github.com/auto-07p/auto-07p
```

Make sure you have FORTRAN installed in your system before running the next step. You can install FORTRAN (in Ubuntu) by typing:

```bash
sudo apt-get install gfortran
```

Then move into the cloned folder and type:

```bash
./configure 
make 
```

Locate the auto.env.sh file. One way to do this in Ubuntu is by typing:

```bash
cd auto-07p
cd cmds
pwd 
```

Append `auto.env.sh` to the directory. Now add the following to the `.bashrc` file:

```bash
source [path-to-auto-07p]/auto.env.sh
```

Make sure that the `AUTO_DIR` variable inside the `auto.env.sh` script point to the same directory you are adding to the `.bashrc` file.

Finally, install auto using pip in order to be able to import AUTO inside the python scripts (If you are using a virtual environment, such as the one provided in this repository, make sure you activate that environment first). Go back to the `auto-07p` folder and type:

```bash
pip install -e .
```

With this you can run `import auto` on the scripts without adding the AUTO path manually.

