function ph_cmb = sense(img,par)
%SENSE combination (for phase).
%   PH_CMB = SENSE(IMG,PAR) combines phase from multiple receivers
%
%   IMG:    raw complex images from multiple receivers, [np nv nv2 ne nrcvrs]
%   PAR:    parameters of the sequence
%   PH_CMB: phase after combination


[~,~,~,ne,nrcvrs] = size(img);
res = par.res;
TE  = par.te + (0:par.ne-1)*par.esp;
TE1 = TE(1);
TE2 = TE(2);


img_diff = img(:,:,:,2,:)./img(:,:,:,1,:);
ph_diff = img_diff./abs(img_diff);
ph_diff_cmb = sum(abs(img(:,:,:,1,:)).*ph_diff,5);

nii = make_nii(angle(ph_diff_cmb),res);
save_nii(nii,'ph_diff.nii');

% perform unwrapping
unix('prelude -p ph_diff.nii -a BET.nii -u unph_diff -m BET_mask.nii -n 8');
unix('gunzip -f unph_diff.nii.gz');

nii = load_nii('unph_diff.nii');
unph_diff_cmb = double(nii.img);

unph_te1_cmb = unph_diff_cmb*TE1/(TE2-TE1);
offsets = img(:,:,:,1,:)./repmat(exp(1j*unph_te1_cmb),[1,1,1,1,4]);
offsets = offsets./abs(offsets);


for chan = 1:nrcvrs
    offsets(:,:,:,:,chan) = smooth3(offsets(:,:,:,:,chan),'box',round(6./res/2)*2+1); 
    offsets(:,:,:,:,chan) = offsets(:,:,:,:,chan)./abs(offsets(:,:,:,:,chan));
end

% combine phase according to complex summation
offsets = repmat(offsets,[1,1,1,ne,1]);
img = img./offsets;
ph_cmb = angle(sum(img,5));
