function [img,params] = swi15_recon(rawfile)
% SWI GRAPPA on 1.5T


% read in the undersampled k-space (0 filled)
% rawfile as a cell of path and filename
% mrdata = siem_read(rawfile);
mrdata = siem_read(rawfile);
k = mrdata.raw;
params = mrdata.params;
sKSpace = params.protocol_header.sKSpace;

% GRAPPA to refill k-space
% Set regularization constant for determination of GRAPPA kernal.
opt.kernreg = 10^-2;

% Set kernal used for GRAPPA (-1 signifies kernal center; 0's are at
% location of skipped lines of k-space)
opt.kern = [1,1,1,1,1;...
    0,0,0,0,0;...
    1,1,1,1,1;...
    0,0,-1,0,0;...
    1,1,1,1,1;...
    0,0,0,0,0;...
    1,1,1,1,1];

% Advanced options
opt.PI_multw = 0;
opt.numwork = 1;

if isfield(params.protocol_header.sPat,'lRefLinesPE')
	% if GARAPPA
	% Extract and reshape reference data
	nref = params.protocol_header.sPat.lRefLinesPE;
	sz = size(k);
	strt = ceil((sz(1)-nref)/2);

	k_ref = k(strt:strt+nref,:,:,:);
	sz_ref = size(k_ref);
	k_ref = reshape(k_ref,[sz_ref(1:3) 1 sz_ref(4)]);
	% size of k_ref: Npe x Nro x NS x (# calibration scans) x RCVRS

	% Calculate weights
	opt.w = grappa_findw(k_ref,opt);
	matlabpool close

	% Reshape k and only take acquired lines
	k_f = permute(k,[1 2 4 3]);
	clear k
	% size of k_f: Npe/2 x Nro x RCVRS x NS x arraydim
	k_f = k_f(1:2:end,:,:,:);

	% Fill lines
	k_f = grappa_fill(k_f,opt);
	matlabpool close

	% Re-fill fully sampled data
	k_f(strt:strt+nref,:,:,:) = permute(k_ref,[1 2 5 3 4]);
	clear k_ref
	% reshape to the final k-space
	k = permute(k_f,[1 2 4 3]);

	% get rid of the last PE line from Corey's recon
	k = k(1:end-1,:,:,:);
end

% partial fourier
sz = size(k);
pf = nan;
flag = sKSpace.ucPhasePartialFourier;
flag(isspace(flag))=[];
switch (flag)
	case {'0x1', 1}  % PF_HALF
		disp('Half fourier is unhandled :(')
	case {'0x2', 2}  % PF_5_8
		pf = 5/8;
	case {'0x4', 4}  % PF_6_8
		pf = 6/8;
	case {'0x8', 8}  % PF_7_8
		pf = 7/8;
	case {'0x10', 10} % PF_OFF
	case {'0x20', 20} % PF_AUTO
	case {'0x16', 16} % "none" in VD??
end
if ~isnan(pf)
	disp(sprintf('Partial Fourier: %d/8', pf*8));
	k = padarray(k, round(sz(1)*(1/(pf)-1)), 'pre');
end

% POCS
if (size(k,4) > 1) && (~isnan(pf)) 
	[~, kspFull] = pocs(permute(k,[4 1 2 3]),20);
	k = permute(kspFull,[2 3 4 1]);
end


% phase resolution
if (sKSpace.dPhaseResolution ~= 1)
	disp(sprintf('Phase Resolution: %1.2f', sKSpace.dPhaseResolution));
	sz = size(k);
	pad_size = round(sz(1)*(1/sKSpace.dPhaseResolution-1));
	k = padarray(padarray(k,round(pad_size/2),'post'),pad_size-round(pad_size/2),'pre');
	% k = padarray(k, round((sz(1)-1)*(1/sKSpace.dPhaseResolution-1)/2));
end


% Asymmetric echo 
sz = size(k);
pad = sKSpace.lBaseResolution - sz(2)/2;
if pad
	disp(sprintf('Asymmetric echo: adding %d lines', pad*2));
	k = padarray(k, [0 pad*2], 'pre');
end


% slice resolution
if (sKSpace.dSliceResolution ~= 1)
	disp(sprintf('Slice Resolution: %1.2f', sKSpace.dSliceResolution));
	sz = size(k);
	pad_size = round(sz(3)*(1/sKSpace.dSliceResolution-1));
	k = padarray(padarray(k,[0,0,round(pad_size/2)],'post'), [0,0,pad_size-round(pad_size/2)],'pre');
	% k = padarray(k, [0, 0, round(sz(3)*(1/sKSpace.dSliceResolution-1)/2)]);
end



% %% remove oversampling in k-space
% k = ifft(k,[],2);
% k(:,size(k,2)/4+1:size(k,2)/4*3,:,:) = [];
% k = fft(k,[],2);

% FFT to image space
img = fftshift(fftshift(fftshift ...
	(fft(fft(fft(fftshift(fftshift(fftshift ...
		(k,1),2),3),[],1),[],2),[],3),1),2),3);

% remove oversampling in RO dimension (2nd dimension)
img = img(:,size(img,2)/4+1:size(img,2)/4*3,:,:);
