function ph_cmb = sense_nof(img,par)
%SENSE combination (for phase). without filter smoothing the offsets
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


% (1) unwrap offset for each channel
for chan = 1:nrcvrs
    nii = make_nii(angle(offsets(:,:,:,1,chan)),res);
    save_nii(nii,['offset_ch' num2str(chan) '.nii']);
    nii = make_nii(abs(img(:,:,:,1,chan)),res);
    save_nii(nii,['magTE1_ch' num2str(chan) '.nii']);
end

bash_command = sprintf(['for ph in offset_ch*\n' ...
'do\n' ...
'       base=`basename $ph`;\n' ...
'       mag="magTE1_ch"${base:9};\n' ...
'       unph="unoffset_ch"${base:9};\n' ...
'       prelude -a $mag -p $ph -u $unph -m BET_mask.nii -n 8&\n' ...
'done\n' ...
'wait\n' ...
'gunzip -f unoffset_ch*.gz\n']);

unix(bash_command);

unoffsets = zeros(size(offsets));
med_unoffsets = unoffsets;
for chan = 1:nrcvrs
    nii = load_nii(['unoffset_ch' num2str(chan) '.nii']);
    unoffsets(:,:,:,1,chan) = double(nii.img);
    med_unoffsets(:,:,:,1,chan) = medfilt3(squeeze(unoffsets(:,:,:,1,chan)),4./res+1);
    nii = make_nii(med_unoffsets(:,:,:,1,chan),res);
    save_nii(nii,['medfilt3_offsets_ch' num2str(chan) '.nii']);
end

% (2) smooth in complex format
%for chan = 1:nrcvrs
%    offsets(:,:,:,:,chan) = smooth3(offsets(:,:,:,:,chan),'box',round(5./res/2)*2+1); % 5mm size kernel
%    offsets(:,:,:,:,chan) = offsets(:,:,:,:,chan)./abs(offsets(:,:,:,:,chan));
%end

% combine phase according to complex summation
med_unoffsets = repmat(exp(1j*med_unoffsets),[1,1,1,ne,1]);
img = img./med_unoffsets;
ph_cmb = angle(sum(img,5));
