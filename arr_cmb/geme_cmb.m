function ph_cmb = geme_cmb(img,vox,te)
%Gradient-echo multi-echo combination (for phase).
%   PH_CMB = GEME_CMB(IMG,PAR) combines phase from multiple receivers
%
%   IMG:    raw complex images from multiple receivers, 5D: [3D_image, echoes, receivers]
%   TE :    echo times
%   vox:    spatial resolution/voxel size, e.g. [1 1 1] for isotropic
%   PH_CMB: phase after combination


[~,~,~,ne,nrcvrs] = size(img);
TE1 = te(1);
TE2 = te(2);

img_diff = img(:,:,:,2,:)./img(:,:,:,1,:);
ph_diff = img_diff./abs(img_diff);
ph_diff_cmb = sum(abs(img(:,:,:,1,:)).*ph_diff,5);
ph_diff_cmb(isnan(ph_diff_cmb)) = 0;

nii = make_nii(angle(ph_diff_cmb),vox);
save_nii(nii,'ph_diff.nii');

% perform unwrapping
unix('prelude -p ph_diff.nii -a BET.nii -u unph_diff -m BET_mask.nii -n 12');
unix('gunzip -f unph_diff.nii.gz');

nii = load_nii('unph_diff.nii');
unph_diff_cmb = double(nii.img);

unph_te1_cmb = unph_diff_cmb*TE1/(TE2-TE1);
offsets = img(:,:,:,1,:)./repmat(exp(1j*unph_te1_cmb),[1,1,1,1,nrcvrs]);
offsets = offsets./abs(offsets);
offsets(isnan(offsets)) = 0;

for chan = 1:nrcvrs
    offsets(:,:,:,:,chan) = smooth3(offsets(:,:,:,:,chan),'box',round(10./vox/2)*2+1); 
    offsets(:,:,:,:,chan) = offsets(:,:,:,:,chan)./abs(offsets(:,:,:,:,chan));
end

%nii = make_nii(angle(offsets),vox);
%save_nii(nii,'offsets.nii');

% combine phase according to complex summation
offsets = repmat(offsets,[1,1,1,ne,1]);
img = img./offsets;
ph_cmb = angle(sum(img,5));

ph_cmb(isnan(ph_cmb)) = 0;
