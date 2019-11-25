% readme

% RECON SUBJECT2 AS AN EXAMPLE


% if have calibration scan p-file

% (1) recon p-file and combine coils using ASSET calibration scan
pfilePath = '/Users/uqhsun8/DATA/MEBRAVO/HONGFU_PRODUCTION_1/Subject2/p-files/P54272.7';
calibrationPfile = '/Users/uqhsun8/DATA/MEBRAVO/HONGFU_PRODUCTION_1/Subject2/p-files/P53248.7';
% kacq_file = 'kacq_yz.txt.1024150629';
outputDir = '/Users/uqhsun8/DATA/MEBRAVO/HONGFU_PRODUCTION_1/Subject2/recon_combined';
% recon_arc_asset(pfilePath, calibrationPfile, kacq_file, outputDir);
recon_arc_asset(pfilePath, calibrationPfile, [], outputDir);

% (2) generate QSM
realDicomsDir = [outputDir '/DICOMs_real'];
imagDicomsDir = [outputDir '/DICOMs_imag'];
path_out = outputDir;
QSM_SPGR_GE_Dir = '/Users/uqhsun8/DATA/MEBRAVO/HONGFU_PRODUCTION_1/Subject2/J_B/Devel_Hba_17109/QSM_SPGR_GE';
qsm_bravo(realDicomsDir, imagDicomsDir, path_out, [], QSM_SPGR_GE_Dir);
% if no SPGR scan then
% qsm_bravo(realDicomsDir, imagDicomsDir, path_out);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if have NO calibration scan p-file or prefer to not use ASSET recon

% (1) recon p-file and output uncombine coils
pfilePath = '/Users/uqhsun8/DATA/MEBRAVO/HONGFU_PRODUCTION_1/Subject2/p-files/P54272.7';
% kacq_file = 'kacq_yz.txt.1024150629';
outputDir = '/Users/uqhsun8/DATA/MEBRAVO/HONGFU_PRODUCTION_1/Subject2/recon_uncombined';
recon_arc(pfilePath, [], outputDir);

% (2) generate QSM
realDicomsDir = [outputDir '/DICOMs_real_uncombined'];
imagDicomsDir = [outputDir '/DICOMs_imag_uncombined'];
path_out = outputDir;
QSM_SPGR_GE_Dir = '/Users/uqhsun8/DATA/MEBRAVO/HONGFU_PRODUCTION_1/Subject2/J_B/Devel_Hba_17109/QSM_SPGR_GE';
qsm_bravo_uncombined(realDicomsDir, imagDicomsDir, path_out, [], QSM_SPGR_GE_Dir);
% if no SPGR scan then
% qsm_bravo_uncombined(realDicomsDir, imagDicomsDir, path_out);


