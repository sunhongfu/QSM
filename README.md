# Welcome to QSM
The repository is for reconstructing Quantitative Susceptiblity Mapping (QSM) images from MRI. The codes cover image recons for single-echo or multiple-echo GRE sequence as well as gradient EPI sequence.

## Recon flow
1. Extract complex img 
2. combine coils 
3. extract brain mask 
4. unwrap phase maps 
5. linearly fit unwrapped phase with TE if multiple echo 
6. background field removal 
7. dipole inversion
	
## Release Note
1. QSM 1.0 (2015-11-23)
  * Recon codes for 1.5T, 3T and 4.7T QSM, including sequences such as SWI, R2s and EPI.
	

## Manual
* Codes location:
  + stable (master) branch: 129.128.117.89\hongfu\Documents\MATLAB\qsm_stable
  + testing (develop) brach: 129.128.117.89\hongfu\Documents\MATLAB\qsm_testing

* Directory structure
  - *15*: recon codes for 1.5T sequences, e.g. EPI (fMRI) and SWI
    + **qsm_epi15.m**
    + **qsm_swi15.m**
  - *47*: recon codes for 4.7T sequences, e.g. EPI (fMRI), SWI and R2*
    + **qsm_epi47.m**
    + **qsm_swi47.m**
    + **qsm_r2s47.m**
  - *PRISMA*: recon codes for 3T PRISMA sequences, e.g. EPI (fMRI), SWI and R2*
    + **qsm_epi_prisma.m**
    + **qsm_swi_prisma.m**
    + **qsm_r2s_prisma.m**
  - *GE*: recon codes for 3T GE sequence, e.g. R2*
    + **qsm_spgr_ge.m**
  - *arr_cmb*: coils combination related codes
    + **adaptive_cmb.m**: adaptive filter method for single-echo, e.g. EPI, SWI
    + **geme_cmb.m**: dual echo approach for multi-echo, e.g. R2*
  - *bkg_rm*: background field removal, including RESHARP/SHARP/ESHARP/PDF/LBV
    + **sharp.m**: SHARP
    + **resharp.m**: RESHARP
    + **projectionontodipolefields.m**: PDF
    + **extendharmonicfield.m**: ESHARP
    + **LBV.m**: LBV
    + **poly2d.m**: 2nd order 2D polynomial fit
    + **poly3d.m**: 2nd order 3D polynomial fit
  - *dip_inv*: dipole inversion with TV regularization
    + **tvdi.m**: total variation dipole inversion
  - *Misc*: other functions including NIFTI and Ryan's small functions

* Usage
  - Call the main QSM function corresponding to the sequence, e.g. `qsm_r2s47` is the function for QSM recon of R2* at 4.7T.
  - Function inputs are 
    + Directory of the raw data for 1.5T/4.7T or directories of both magnitude and unfiltered phase DICOMs for PRISMA/GE 3T)
    + User defined directory for QSM output results
    + User specified parameters "options"
  - Examples:
  
        ```Matlab
        options.bkg_rm='resharp';
        options.ph_unwrap='laplacian';
        qsm_swi47('FID_DIR','OUTPUT_DIR',options);
        ```
        
  - For other advanced usage, see help, e.g. `help qsm_swi_prisma`

* Common user-changed options:
  - **bet_thr**: threshold level for BET extracting the brain mask, by default is *0.3-0.5* depending on the sequence, smaller threshold keeps more region of the brain
  - **ph_unwrap**: phase unwrapping methods, can be *'prelude'*, *'laplacian'* or *'bestpath'*
  - **bkg_rm**: background field removal methods, can be *'sharp'*, *'pdf'*,*'resharp'*,*'esharp'* or *'lbv'*, can pick multiple methods to compare, e.g. `options.bkg_rm={'resharp','lbv'}`
  - **smv_rad**: radius in mm of SHARP/RESHARP/ESHARP kernel (erosion size, ESHARP recovers some)
  - **tik_reg**: tikhonov regularization for RESHARP, by default *1e-3*, bigger value more regularization
  - **tv_reg**: total variation regularization, by default *5e-4*, bigger value gives smoother result
