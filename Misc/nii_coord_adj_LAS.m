function nii_coord_adj_LAS(nii_input_name, nii_ref_name)

% data is stored in LAS instead of RAS
% this is to be compatible with some nifti that data is stored in LAS, and header displays it back to RAS

nii_input = load_untouch_nii(nii_input_name);
nii_ref = load_untouch_nii(nii_ref_name);

nii_output = nii_ref;
nii_output.img = (flipdim((single(nii_input.img)),2));

nii_output.hdr.dime.datatype = 16;
nii_output.hdr.dime.bitpix = 32;

nii_output.hdr.dime.glmax = max(nii_input.img(:));
nii_output.hdr.dime.glmin = min(nii_input.img(:));
nii_output.hdr.dime.dim = nii_input.hdr.dime.dim;

[filepath,name,ext] = fileparts(nii_input_name);
nii_output_name = fullfile(filepath,[name '_adj'],ext);

save_untouch_nii(nii_output,nii_output_name);


