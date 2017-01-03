# Welcome to QSM
The repository is for reconstructing Quantitative Susceptiblity Mapping (QSM) images from MRI. The codes cover image recons for single-echo SWI or multiple-echo GRE sequence as well as gradient EPI sequence.

## Recon flow
1. Extract complex img 
2. combine coils 
3. extract brain mask 
4. unwrap phase maps 
5. linearly fit unwrapped phase with TE if multiple echo 
6. background field removal 
7. dipole inversion
	
## Release Note
1. QSM 2.0 (2017-01-01)
  * Recon codes for 1.5T, 3T and 4.7T QSM, including sequences such as SWI, R2* and EPI.
	

## Manual
* Directory structure
  - *Siemens_1p5T*: recon codes for Siemens 1.5T sequences, e.g. EPI (fMRI) and SWI
    + **qsm_epi15.m**
    + **qsm_swi15.m**
    + **qsm_hemo15.m**
  - *Varian_4p7T*: recon codes for Varian 4.7T sequences, e.g. EPI (fMRI), SWI and R2*
    + **qsm_epi47.m**
    + **qsm_swi47.m**
    + **qsm_r2s47.m**
  - *Siemens_3T*: recon codes for 3T PRISMA sequences, e.g. EPI (fMRI), SWI and R2*
    + **qsm_epi_prisma.m**
    + **qsm_swi_prisma.m**
    + **qsm_r2s_prisma.m**
    + **qsm_hemo_prisma.m**
  - *GE_3T*: recon codes for 3T GE sequence, e.g. R2*
    + **qsm_spgr_ge.m**
  - *coil_combination*: coils combination related codes
    + **adaptive_cmb.m**: adaptive filter method for single-echo, e.g. EPI, SWI
    + **geme_cmb.m**: dual echo approach for multi-echo, e.g. R2*
  - *background_field_removal*: background field removal, including RESHARP/SHARP/ESHARP/PDF/LBV
    + **sharp.m**: SHARP
    + **resharp.m**: RESHARP
    + **projectionontodipolefields.m**: PDF
    + **extendharmonicfield.m**: ESHARP
    + **LBV.m**: LBV
    + **poly2d.m**: 2nd order 2D polynomial fit
    + **poly3d.m**: 2nd order 3D polynomial fit
  - *dipole_inversion*: dipole inversion with TV regularization
    + **tvdi.m**: total variation dipole inversion
    + **tikhonov_qsm.m**: Tikhonov total field inversion
  - *Misc*: other functions including NIFTI and other small functions

* Usage
  - Call the main QSM function corresponding to the sequence, e.g. `qsm_r2s47` is the function for QSM recon of R2* at 4.7T.
  - Function inputs are 
    + Directory of the raw data for 1.5T/4.7T or directories of both magnitude and unfiltered phase DICOMs for Siemens/GE 3T)
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
