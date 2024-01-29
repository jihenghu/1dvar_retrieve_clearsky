# rttov_fwd_k_example

To run the RTTOV in built Kmatrix test, offine.

## build your own
1. replace all the RTTOV-related path to $your_RTTOV_Intall_PATH in `Makefile.sh` ; 
2. specify your own `$COEF_FILENAME` in `run_example_k.sh`:
```bash
COEF_FILENAME=$PATH_TO_Your_RTTOV"/rttov13pred54L/rtcoef_gpm_1_gmi.dat"
```

## Compile  
```bash
> sh Makefile.sh
```

## Run 
```bash
> sh run_example_k.sh
```