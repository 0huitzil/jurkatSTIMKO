# A model for the role of STIM1 and STIM2 in Ca2+ oscillations in T-lymphocytes

## Installing AUTO 
First, clone the AUTO repository. Open a terminal in the folder you wish to put the AUTO files in and type:
```
git clone https://github.com/auto-07p/auto-07p
```
Make sure you have FORTRAN installed in your system before running the next step. You can install fortran (in Ubuntu) by typing
```
sudo apt-get install gfortran
```
Then move into the cloned folder and type:
```
./configure 
make 
```
Locate the auto.env.sh file. One way to do this in Ubuntu is by typing
```
cd auto-07p
cd cmds
pwd 
```
and append `auto.env.sh` to the directory. Now add the following to the `.bashrc` file
```
source [path-to-auto-07p]/auto.env.sh
```
Make sure that the AUTO_DIR variable inside the `auto.env.sh` script point to the same directory you are adding to the `.bashrc` file. 

Finally, in order to be able to import AUTO inside the python scripts, go back to the `auto-07p` folder and type 
```
pip install -e .
```
With this you can run `import auto` on the scripts without adding the AUTO path manually. 