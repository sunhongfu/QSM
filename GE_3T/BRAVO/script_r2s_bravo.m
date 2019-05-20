% MATLAB script for computing R2* map from multi-echo BRAVO sequence
load('all.mat','mag','TE','mask','imsize','vox');
[R2 T2 amp] = r2imgfit2(double(mag),TE,repmat(mask,[1 1 1 imsize(4)]));
nii = make_nii(R2,vox);
save_nii(nii,'R2.nii');
nii = make_nii(T2,vox);
save_nii(nii,'T2.nii');
nii = make_nii(amp,vox);
save_nii(nii,'amp.nii');

