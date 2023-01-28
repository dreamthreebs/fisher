# Fishmat
Fisher matrix method to calculate constraint on dark matter parameters
## Dependency
ifort(intel fortran) or gfortran, icc(intel c++), python(numpy, matplotlib, scipy)
## Structure
./camb/calc/ : all codes for calculation partial derivative(pd), and one error sigma(sig)  
./camb/data/ : pd data calculated by getpd.py and sigma data calculated by getsig.py  
./camb/input/ : noise level(nls) and foreground residual(fgres) for ali, pico, and planck  
./camb/plot/ : plot all figures and save in figure folder  
## How to use
make Hyrec; make camb  
run ./camb/calc/getpd.py to get all pd data in ./camb/data/pd\_data/\*CLprime.npy  
run ./camb/calc/getsig.py to get all one sigma error in ./camb/data/sigma\_data/sigma\*.npy  
run ./camb/plot/plot\*.py to get different figures of what you want  

Precautions: 1.intel: remember to set
```bash
source /path/to/intel/oneapi/setvars.sh
source /path/to/intel/oneapi/compiler/your intel version/env/vars.sh
```
2.before using codes or do different research on scattering cross section,you need to
```bash
cd Hyrec;make clean;rm libhyrec.a;make
cd camb;make clean;rm camb;make
./camb your params.ini
```
3.sometimes noise and fgres varies more frequently, so the pd can be calculated first. Then carry on calculations on sigma
