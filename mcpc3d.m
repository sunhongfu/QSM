function ph_cmb = mcpc3d(img,par)
%MCPC-3D combination (for phase).
%   IMG_CMB = mcpc3d(IMG,PAR) combines data from multiple receivers
%
%   IMG: raw complex images from multiple receivers, [np nv nv2 ne nrcvrs]
%   PAR: parameters of the sequence
%
%   requires: bash script 'unwrap'


[np,nv,nv2,ne,nrcvrs] = size(img);
res = par.res;
TE = par.te + (0:par.ne-1)*par.esp;
TE1 = TE(1);
TE2 = TE(2);

% unwrap first 2 echoes, all channels
for echo = 1:2
    for chan = 1:nrcvrs
        nii = make_nii(abs(img(:,:,:,echo,chan)));
        save_nii(nii,['mag_te' num2str(echo) '_ch' num2str(chan) '.nii']);
	nii = make_nii(angle(img(:,:,:,echo,chan)));
	save_nii(nii,['ph_te' num2str(echo) '_ch' num2str(chan) '.nii']);
    end
end

% perform unwrapping (outputs: unph_te*_chan*.nii)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unix('~/projects/QSM/scripts/unwrap ph_te*_ch*.nii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unph = zeros([np,nv,nv2,2,nrcvrs]);
for echo = 1:2
    for chan = 1:nrcvrs
	nii = load_nii(['unph_te' num2str(echo) '_ch' num2str(chan) '.nii']);
	unph(:,:,:,echo,chan) = double(nii.img);
    end
end

% calculate phase-offsets according to dual-echo approach
offsets = squeeze((TE1*unph(:,:,:,2,:)-TE2*unph(:,:,:,1,:))/(TE1-TE2));

% smooth phase-offsets using median filter 5*5*5 pixels
for chan = 1:nrcvrs
    nii = make_nii(offsets(:,:,:,chan));
    save_nii(nii,['raw_offset' num2str(chan) '.nii']);
    offsets(:,:,:,chan) = medfilt3(offsets(:,:,:,chan), [5,5,5]);
    nii = make_nii(offsets(:,:,:,chan));
    save_nii(nii,['smooth_offset' num2str(chan) '.nii']);
end

% combine phase according to complex summation
offsets = permute(repmat(offsets,[1,1,1,1,ne]),[1,2,3,5,4]);
img = img./exp(1j*offsets);
ph_cmb = angle(sum(img,5));
