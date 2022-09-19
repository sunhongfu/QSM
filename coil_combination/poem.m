function [ph_cmb,mag_cmb] = poem(mag, pha, vox, te, mask, smooth_method, parpool_flag)
%Gradient-echo multi-echo combination (for phase).
%   PH_CMB = POEM(MAG, PHA, VOX, TE, MASK, SMOOTH_METHOD) combines phase from multiple receivers
%
%   MAG/PHA:     raw complex images from multiple receivers, 5D: [3D_image, echoes, receiver channels]
%   TE :     echo times
%   MASK:    brain mask
%   VOX:     spatial resolution/voxel size, e.g. [1 1 1] for isotropic
%   PH_CMB:  phase after combination
%   MAG_CMB: magnitude after combination
%   SMOOTH:  smooth methods(1) smooth3, (2) poly3, (3) poly3_nlcg, (4) gaussian

if ~ exist('smooth_method','var') || isempty(smooth_method)
    smooth_method = 'smooth3';
end

if ~ exist('parpool_flag','var') || isempty(parpool_flag)
    parpool_flag = 1;
end

if isdeployed
    parpool_flag = 0;
end

[~,~,~,ne,nrcvrs] = size(mag);
TE1 = te(1);
TE2 = te(2);
imsize = size(mag);

ph_diff = exp(1j*pha(:,:,:,2,:))./exp(1j*pha(:,:,:,1,:)) ;
ph_diff_cmb = sum(ph_diff,5);
ph_diff_cmb(isnan(ph_diff_cmb)) = 0;

nii = make_nii(angle(ph_diff_cmb),vox);
save_nii(nii,'ph_diff.nii');

% % perform unwrapping
% method (1)
% unix('prelude -p ph_diff.nii -a BET.nii -u unph_diff -m BET_mask.nii -n 12');
% unix('gunzip -f unph_diff.nii.gz');
% nii = load_nii('unph_diff.nii');
% unph_diff_cmb = double(nii.img);
mag1 = sqrt(sum((mag(:,:,:,1,:).^2),5));
mask_input = mask;
mask = (mag1 > 0.1*median(mag1(logical(mask(:)))));
mask = mask | mask_input;

% method (2)
% best path unwrapping
[pathstr, ~, ~] = fileparts(which('3DSRNCP.m'));
setenv('pathstr',pathstr);
setenv('nv',num2str(imsize(1)));
setenv('np',num2str(imsize(2)));
setenv('ns',num2str(imsize(3)));

fid = fopen('wrapped_phase_diff.dat','w');
fwrite(fid,angle(ph_diff_cmb),'float');
fclose(fid);

mask_unwrp = uint8(mask*255);
fid = fopen('mask_unwrp.dat','w');
fwrite(fid,mask_unwrp,'uchar');
fclose(fid);

if isdeployed
    bash_script = ['~/bin/3DSRNCP wrapped_phase_diff.dat mask_unwrp.dat ' ...
    'unwrapped_phase_diff.dat $nv $np $ns reliability_diff.dat'];
else
    bash_script = ['${pathstr}/3DSRNCP wrapped_phase_diff.dat mask_unwrp.dat ' ...
    'unwrapped_phase_diff.dat $nv $np $ns reliability_diff.dat'];
end
unix(bash_script) ;

fid = fopen('unwrapped_phase_diff.dat','r');
tmp = fread(fid,'float');
unph_diff_cmb = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
fclose(fid);

nii = make_nii(unph_diff_cmb,vox);
save_nii(nii,'unph_diff.nii');


% calculate initial phase offsets
unph_te1_cmb = unph_diff_cmb*TE1/(TE2-TE1);
offsets = exp(1j*pha(:,:,:,1,:))./repmat(exp(1j*unph_te1_cmb),[1,1,1,1,nrcvrs]);
offsets(isnan(offsets)) = 0;

nii = make_nii(angle(offsets),vox);
save_nii(nii,'offsets_raw.nii');
% 
% smooth offsets
% maybe later change to smooth the real and imag parts seperately, and try
% guassian filter!
if strcmpi('smooth3',smooth_method)
    if parpool_flag
        parpool;
        parfor chan = 1:nrcvrs
            offsets(:,:,:,1,chan) = smooth3(offsets(:,:,:,1,chan),'box',round(5)*2+1); 
    %       offsets(:,:,:,1,chan) = smooth3(offsets(:,:,:,1,chan),'box',round(2./vox)*2+1); 
            offsets(:,:,:,1,chan) = offsets(:,:,:,1,chan)./abs(offsets(:,:,:,1,chan));
        end
        delete(gcp('nocreate'));
    else
        for chan = 1:nrcvrs
            offsets(:,:,:,1,chan) = smooth3(offsets(:,:,:,1,chan),'box',round(5)*2+1); 
    %       offsets(:,:,:,1,chan) = smooth3(offsets(:,:,:,1,chan),'box',round(2./vox)*2+1); 
            offsets(:,:,:,1,chan) = offsets(:,:,:,1,chan)./abs(offsets(:,:,:,1,chan));
        end
    end

elseif strcmpi('gaussian',smooth_method)
    if parpool_flag
        parpool;
        parfor chan = 1:nrcvrs
            offsets(:,:,:,1,chan) = imgaussfilt3(real(offsets(:,:,:,1,chan)),4) + 1j*imgaussfilt3(imag(offsets(:,:,:,1,chan)),4);
            offsets(:,:,:,1,chan) = offsets(:,:,:,1,chan)./abs(offsets(:,:,:,1,chan));
        end
        delete(gcp('nocreate'));
    else
        for chan = 1:nrcvrs
            offsets(:,:,:,1,chan) = imgaussfilt3(real(offsets(:,:,:,1,chan)),4) + 1j*imgaussfilt3(imag(offsets(:,:,:,1,chan)),4);
            offsets(:,:,:,1,chan) = offsets(:,:,:,1,chan)./abs(offsets(:,:,:,1,chan));
        end
    end

elseif strcmpi('poly3',smooth_method)
    for chan = 1:nrcvrs
        fid = fopen(['wrapped_offsets_chan' num2str(chan) '.dat'],'w');
        fwrite(fid,angle(offsets(:,:,:,chan)),'float');
        fclose(fid);
        setenv('chan',num2str(chan));
        if isdeployed
            bash_script = ['~/bin/3DSRNCP wrapped_offsets_chan${chan}.dat mask_unwrp.dat ' ...
            'unwrapped_offsets_chan${chan}.dat $nv $np $ns reliability_diff.dat'];
        else
            bash_script = ['${pathstr}/3DSRNCP wrapped_offsets_chan${chan}.dat mask_unwrp.dat ' ...
            'unwrapped_offsets_chan${chan}.dat $nv $np $ns reliability_diff.dat'];
        end
        unix(bash_script) ;
        fid = fopen(['unwrapped_offsets_chan' num2str(chan) '.dat'],'r');
        tmp = fread(fid,'float');
        unph_offsets(:,:,:,chan) = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
        fclose(fid);
        offsets(:,:,:,chan) = poly3d(unph_offsets(:,:,:,chan),mask,3);
    end
    offsets = exp(1j*offsets);
elseif strcmpi('poly3_nlcg',smooth_method)
    for chan = 1:nrcvrs
        offsets(:,:,:,:,chan) = poly3d_nonlinear(offsets(:,:,:,:,chan),mask,3);
    end
else
    error('what method to use for smoothing? smooth3 or poly3 or poly3_nlcg')
end
nii = make_nii(angle(offsets),vox);
save_nii(nii,'offsets_smooth.nii');


% combine phase according to complex summation
offsets = repmat(offsets,[1,1,1,ne,1]);
img_cmb = sum(mag.*exp(1j*pha)./offsets,5);
img_cmb(isnan(img_cmb)) = 0;
ph_cmb = angle(img_cmb);
ph_cmb(isnan(ph_cmb)) = 0;
mag_cmb = abs(img_cmb);
mag_cmb(isnan(mag_cmb)) = 0;

