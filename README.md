# CCFpams: Atmospheric Stellar Parameters from Cross-Correlation Functions

**`CCFpams` v1.0 by Luca Malavolta - 2017**    

CCFpams is a novel approach that allows the measurement of stellar temperature,
metallicity and gravity within a few seconds and in a completely automated fashion. Rather than performing comparisons with spectral libraries, our technique is based on the determination of several cross-correlation functions (CCFs) obtained by including spectral features with different sensitivity to the photospheric parameters. We use literature stellar parameters of high signal-to-noise (SNR), high-resolution HARPS spectra of FGK Main Sequence stars to calibrate the stellar parameters as a function of CCF areas.  For FGK
stars we achieve a precision of 50K in temperature, 0.09 dex in gravity and 0.035 dex in metallicity at SNR=50 while the precision for observation with SNR>100 and the overall accuracy are constrained by the literature values used to calibrate the CCFs.

The paper describing the technique has been accepted by MNRAS and it will be available on Monday 05/08/2017. In this repository the code to extract the parameters from HARPS and HARPS-N is provided.

## CCFpams user guide

1. [Install and Compile](### Install and Compile)
2. [Prepare the observations](### Prepare the observations)

### Install and Compile
 `CCFpams` is avaialble at [https://github.com/LucaMalavolta/PyORBIT](https://github.com/LucaMalavolta/CCFpams)

 Download the .zip files or clone the repository. `harps_input2pams.f90` and `harpn_input2pams.f90` are the two FORTRAN90 codes to compute the stellar parameters from HARPS and HARPS-N data respectively.

 The subroutines required by the two programs are inside the `routines_f90` folder. Calibration files required by the two programs are stored in the `mask_calib` folder.

 Before compiling, you have to declare the path of the code in the code itself. In this way the code will be able to find automatically all the required calibration files. Additionally you can specify the directory of the HARPS/HARPS-N archive of the _DRS reduced_ files.

 To do so, open `harps_input2pams_v3.f90`/`harpn_input2pams_v3.f90` with a text editor and change `code_path` and (optionally) `archive_harpn` with the full path of your folders.

 ```fortran
  character (len=*), parameter :: &
     code_path = '/home/malavolta/CODE/CCFpams/', &
     archive_harpn = '/home/malavolta/data/HARPN/data/'
  ```

To compile the program, just execute in a shell:

 ```sh
  $ ./COMPILE
 ```
 The script will create two executable files, `harps_input2pams.e` and `harpn_input2pams_v3.e`

### Prepare the observations
