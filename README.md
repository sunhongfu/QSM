# Welcome to QSM
For deep learning-based QSM methods, check out the other github repo: [deepMRI](https://github.com/sunhongfu/deepMRI)

The repository is for reconstructing Quantitative Susceptibility Mapping (QSM) images from MRI. The codes cover image recons for single-echo SWI or multiple-echo GRE sequence as well as gradient EPI sequence.

# References
* For the background field removal method RESHARP, please reference:  
***H. Sun, A.H. Wilman. Background field removal using spherical mean value filtering and Tikhonov regularization. Magn Reson Med. 2014 Mar;71(3):1151-7.***
* For the POEM multi-channel coil combination method, please referece:  
***H. Sun, J.O. Cleary, R. Glarin, S.C. Kolbe, R.J. Ordidge, B.A. Moffat, G.B. Pike; Extracting more for less: Multi-echo MP2RAGE for simultaneous T1-weighted imaging, T1 mapping, R2\* mapping, SWI, and QSM from a single acquisition.***
* For the EPI-QSM processing pipeline, please referece:  
***H. Sun, A.H. Wilman. Quantitative susceptibility mapping using single-shot echo-planar imaging. Magn Reson Med. 2015 May;73(5):1932-8.***
* For the hemorrhage-QSM method and processing pipeline, please referece:  
***H. Sun, M. Kate, L.C. Gioia, D.J. Emery, K. Butcher, A.H. Wilman. Quantitative susceptibility mapping using a superposed dipole inversion method: Application to intracranial hemorrhage. Magn Reson Med. 2016 Sep;76(3):781-91.***

## Recon flow
1. Extract complex img from DICOMs or raw files
2. removal phase-offsets if multiple-echo, then combine coils 
3. extract brain mask 
4. unwrap phase maps 
5. linearly fit unwrapped phase with TE if multiple-echo 
6. background field removal 
7. dipole inversion
	
## Support scanner platforms
  * 1.5 T (Siemens)
  * 3   T (Siemens, GE, and Philips)
  * 4.7 T (Varian)
  * 7   T (Siemens)
  * 9.4 T (Bruker)
  * Support sequences of single-echo SWI, multi-echo GRE and EPI.

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
  - *Philips_3T*: recon codes for 3T Philips sequence, currently only R2*
    + **qsm_spgr_philips.m**
  - *Siemens_7T*: recon codes for 7T Magnetom sequences, e.g. SWI, ME-GRE and ME-MP2RAGE
    + **qsm_swi.m**
    + **mp2rage.m**
    + **qsm_7T_bipolar.m**
    + **qsm_7T_unipolar.m**
  - *Bruker_9p4T*: recon codes for 9.4T Bruker 3D multi-echo sequence
    + **qsm_r2s_3d_bruker.m**  
  - *coil_combination*: coils combination related codes
    + **adaptive_cmb.m**: adaptive filter method for single-echo, e.g. EPI, SWI
    + **poem.m**: POEM (Phase-Offsets Estimation from Multi-echoes) coil combination
  - *background_field_removal*: background field removal, including RESHARP/SHARP/ESHARP/PDF/LBV
    + **sharp.m**: SHARP
    + **resharp.m**: RESHARP method -- can tweak the ker_rad (kernel size) and tik_reg (regularization)
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
    + Directory of the raw data for 1.5T/4.7T or directories of both magnitude and unfiltered phase DICOMs for Siemens/GE/Philips 3T)
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
  - **tik_reg**: tikhonov regularization for RESHARP, by default *1e-4*, bigger value more regularization
  - **tv_reg**: total variation regularization, by default *5e-4*, bigger value gives smoother result
