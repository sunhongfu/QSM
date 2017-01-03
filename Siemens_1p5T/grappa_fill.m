function U = grappa_fill(k,opt,low_int_mask)
%   Fills in missing k-space lines using GRAPPA weightings obtained from
%   grappa_findw.m. This function currently only supports rate 2 GRAPPA.
%   
%   Usage: k = grappa_findw(kr,opt,low_int_mask)
%   Author: Corey Baron
%   Date: 10/2010
%   
%   Input:
%   kr: k-space data for partially sampled data of size
%       Npe/2 x Nro x RCVRS x NS x arraydim
%       This data should already be phase corrected and gridded.
%   opt:    structure with recon options: weights (opt.w), kernal
%           (opt.kern)
%       opt.w: GRAPPA weights of size 
%           (RCVRS * # non-zero kernal entries) x RCVRS x NS x (# calibration
%           scans)
%       opt.kern: kernal used for GRAPPA weights. A value of -1 corresponds to the
%           center of the kernal, 0 for unused points, and 1 for used points. The
%           kernal must include lines that will be missing in the undersampled data
%           (these should be rows of zeros in the kernal). The kernal used here
%           should exactly match the kernal used in grappa_findw.m.
%       opt.numwork: number of workers to use for parallel processing.
%       opt.costtype: only relevant if an optimal set of grappa weights
%            needs to be chosen when more than 1 are supplied. Choices:
%               1 (default) - image entropy
%               2 - total size of a thresholded image
%               3 - total intensity outside the mask obtained from the
%                   fully sampled reference data.
%   low_int_mask: only required if opt.costtype == 3 and more than 1 set 
%       of grappa weights are supplied. This is an image mask created from 
%       the fully sample reference data, so that the function knows where 
%       to look for and suppress ghost signal.
%   
%   Output:
%   k: full k-space data of size
%       Npe x Nro x RCVRS x NS x arraydim

% Check inputs
if nargin < 3
    low_int_mask = [];
else
    low_int_mask = logical(round(1-low_int_mask));
end
if ~isfield(opt,'verbose') || isempty(opt.verbose)
    opt.verbose = 0;
end
verbose = opt.verbose;
if ~isfield(opt,'numwork') || isempty(opt.numwork)
    opt.numwork = 1;
end
ppm1 = [];

% The below flag specifies whether to shift the position of the image to
% find the best match with the reference image. This seems to be required
% mostly for gradient echo EPI.
if ~isfield(opt,'epi_ge') || isempty(opt.epi_ge)
    opt.epi_ge = 0;
end

if isfield(opt,'costtype') && (opt.costtype == 2)
    cstfcn = @costfcn2;
elseif isfield(opt,'costtype') && (opt.costtype == 3)
    cstfcn = @costfcn3;
else
    cstfcn = @costfcn1;
end

hack_outervox = 0;
if hack_outervox==1
    % Note: this is just to check if shortcut I took for outer k-space
    % voxels matters. It doesn't seem to...
    warning('GRAPPA_fill hack enabled')
end

% Define some values
% k = permute(kr,[1 2 5 3 4]);
npe = size(k,1);
nro = size(k,2);
npnts = size(k,1)*size(k,2);
nrcvrs = size(k,3);
nslices = size(k,4);
nad = size(k,5);
classk = class(k);

% Start progress bar
ticID = tic;
if verbose > 0
    ppm1 = ParforProgressStarter2('GRAPPA - finding optimal weights', nslices, 0.1);
end

% Find center of kernal
[Ci, Cj] = find(opt.kern == -1);

% Remove lines of the kernal that correspond to unnaquired lines in
% k-space.
if mod(Ci,2) == 0
    inds = 1:2:size(opt.kern,1);
else
    inds = 2:2:size(opt.kern,1);
end
kern_s = opt.kern(inds,:);
Ci = round(Ci/2);

% Find the non-zero entries in the kernal
[I,J] = find(kern_s>0);
I = I - Ci;
J = J - Cj;
Il = length(I);
mI = max(abs(I));
mJ = max(abs(J));

% Create array to hold index of best weighting for each slice
b_ind = ones(nslices,1);
b_shft = ones(nslices,1);

% Get best weighting factor for each array dimension (use b0 image,
% assuming it's ad=1)
szw4 = size(opt.w,4);
if opt.epi_ge == 1
    % If gradient echo EPI, it's possible reference scan image was shifted
    % wrt undersampled data due to error in frequency offset. So, we look
    % for the best offset. Spans the entire FOV in 1/4 pixel increments.
    % TODO: for distortion match PI reference, I could reduce the span to
    % only a few voxels, and decrease the increment.
    pix = -npe/4:0.25:npe/4;
else
    pix = 0;
end
if (szw4 > 1) || (opt.epi_ge == 1)
    if opt.numwork > 1
        tmp1 = opt.epi_ge;
        tmp2 = opt.w;
        parfor (s = 1:nslices,opt.numwork)
            % Figure out index of best weighting factors, as well as index
            % of best image shift
            [b_ind(s), b_shft(s)] = get_b(k(:,:,:,s,:),tmp1,tmp2(:,:,s,:),pix,Il,I,J,classk,cstfcn,low_int_mask(:,:,s));
            % Update progress bar
            if verbose > 0
                ppm1.increment(s);
            end
        end
        clear tmp1 tmp2
        if opt.epi_ge == 1
            % Use the best shifts to shift the images
            for s=1:nslices
                ph_ramp = exp(sqrt(-1)*pi*pix(b_shft(s))*(-1:2/size(k(:,:,:,s,:),1):1-1/size(k(:,:,:,s,:),1)));
                k(:,:,:,s,:) = k(:,:,:,s,:) .* repmat(ph_ramp.',[1 size(k(:,:,:,s,:),2) size(k(:,:,:,s,:),3) 1 size(k(:,:,:,s,:),5)]);
            end
        end
    else
        for s = 1:nslices
            % Figure out index of best weighting factors, as well as index
            % of best image shift
            [b_ind(s), b_shft(s)] = get_b(k(:,:,:,s,:),opt.epi_ge,opt.w(:,:,s,:),pix,Il,I,J,classk,cstfcn,low_int_mask(:,:,s));
            
            if opt.epi_ge == 1
                ph_ramp = exp(sqrt(-1)*pi*pix(b_shft(s))*(-1:2/size(k(:,:,:,s,:),1):1-1/size(k(:,:,:,s,:),1)));
                k(:,:,:,s,:) = k(:,:,:,s,:) .* repmat(ph_ramp.',[1 size(k(:,:,:,s,:),2) size(k(:,:,:,s,:),3) 1 size(k(:,:,:,s,:),5)]);
            end
            % Update progress bar
            if verbose > 0
                ppm1.increment(s);
            end
        end
    end
end

% Finish up the progress bar
if verbose > 0
    try % use try / catch here, since delete(struct) will raise an error.
      delete(ppm1);
    catch me %#ok<NASGU>
    end
    dur = toc(ticID);
    remove_summary(dur);
    % Start next progress bar
    ppm1 = ParforProgressStarter2('GRAPPA line filling progress', nad*nslices+nad*nslices*nrcvrs+nslices, 0.1);
end

% Find the missing lines using best weighting factors.
% Create array to hold calculated lines.
U = zeros([npnts,nrcvrs,nslices,nad],class(k)) +...
    1i*ones([npnts,nrcvrs,nslices,nad],class(k));
Usz = size(U);
if opt.numwork > 0
    w = opt.w;
    for a=1:nad
        k_a = k(:,:,:,:,a);
        parfor (s = 1:nslices, opt.numwork)
            U_a = findmissing(Usz,k_a(:,:,:,s),I,J,w,b_ind,a,s);
            U(:,:,s,a) = U_a;
            % Update progress bar
            if verbose > 0
                ppm1.increment(a*(nslices-1)+s);
            end
        end
    end
    clear w k_a
else
    for a=1:nad
        for s = 1:nslices
            U_a = findmissing(Usz,k(:,:,:,s,a),I,J,opt.w,b_ind,a,s);
            U(:,:,s,a) = U_a;
            % Update progress bar
            if verbose > 0
                ppm1.increment(a*(nslices-1)+s);
            end
        end
    end
end

% Combine k and U for memory efficiency.
U = reshape(U,[npe nro nrcvrs nslices nad]);
U = [k;U];
clear k

% Interleave the odd and even lines of k-space.
if opt.numwork > 1
    parfor (a=1:nad*nslices*nrcvrs, opt.numwork)
        tmp = U(:,:,a);
        tmp2 = tmp(1:npe,:);
        tmp(2:2:end,:) = tmp((npe+1):end,:);
        tmp(1:2:end,:) = tmp2;
        U(:,:,a) = tmp;
%         U(:,:,a) = upsample(tmp(1:npe,:),2,0) + upsample(tmp((npe+1):end,:),2,1);
        % Update progress bar
        if verbose > 0
            ppm1.increment(a+nad*nslices);
        end
    end
else
    for a=1:nad*nslices*nrcvrs
        tmp = U(:,:,a);
        tmp2 = tmp(1:npe,:);
        tmp(2:2:end,:) = tmp((npe+1):end,:);
        tmp(1:2:end,:) = tmp2;
        U(:,:,a) = tmp;
%         U(:,:,a) = upsample(U(1:npe,:,a),2,0) + upsample(U((npe+1):end,:,a),2,1);
        % Update progress bar
        if verbose > 0
            ppm1.increment(a+nad*nslices);
        end
    end
end
clear tmp*

% Return slice position back to original
if opt.epi_ge == 1
    for s = 1:nslices
        ph_ramp = exp(-sqrt(-1)*pi*pix(b_shft(s))*(-1:2/size(U,1):1-1/size(U,1)));
        U(:,:,:,s,:) = U(:,:,:,s,:) .* repmat(ph_ramp.',[1 size(U,2) size(U,3) 1 size(U,5)]);
        % Update progress bar
        if verbose > 0
            ppm1.increment(s+nad*nslices+nad*nslices*nrcvrs);
        end
    end
end

% Hack: eliminate the outer voxels in k-space that were incorrectly
% calculated.
if hack_outervox == 1
    U([1:mI,end-mI+1:end],[1:mJ,end-mJ+1:end],:,:,:) = 0;
end

% Finish up the progress bar
if verbose > 0
    try % use try / catch here, since delete(struct) will raise an error.
      delete(ppm1);
    catch me %#ok<NASGU>
    end
    dur = toc(ticID);
    remove_summary(dur);
end

end

function U_a = findmissing(Usz,k_a,I,J,w,b_ind,a,s)
% Function to apply the GRAPPA kernal to find missing lines of k-space

Il = length(I);
npnts = size(k_a,1)*size(k_a,2);
nrcvrs = size(k_a,3);
A = zeros([npnts,Il*nrcvrs],class(k_a)) +...
            1i*ones([npnts,Il*nrcvrs],class(k_a));

% Use circular shifts to isolate the desired elements of k
for r = 1:Usz(2)
    for n=1:Il
        A(:,(r-1)*Il+n) =...
            reshape( circshift(k_a(:,:,r),[-I(n),-J(n)]),[],1);
    end
end

% Calculate missing lines using best weighting factors
U_a = A*w(:,:,s,b_ind(s));

end

function [b_ind_sl, b_shft_sl] = get_b(k_sl,epi_ge,w,pix,Il,I,J,classk,cstfcn,low_int_mask)
% Function to figure out index of best weighting factors (b_ind_sl), as 
% well as index of best image shift (b_shft_sl), for a slice.

% Make some definitions and initializations
N = length(pix);
szw4 = size(w,4);
npe = size(k_sl,1);
nro = size(k_sl,2);
nrcvrs = size(k_sl,3);
npnts = npe*nro;
A_t = zeros([npnts,Il*nrcvrs],classk) +...
    1i*ones([npnts,Il*nrcvrs],classk);

% Choose which array dimension to use for getting optimal weighting
% factors. Typically for DTI the first one is the b0, which is the
% best one to use.
ad = 1;

ent = zeros(szw4,N);
all_flag = 0;
% Loop through all image shifts
for n2=1:N
    
    % Apply image shift
    tmp = k_sl;
    tmp = tmp(:,:,:,1,ad);
    if epi_ge == 1
        ph_ramp_t = exp(sqrt(-1)*pi*pix(n2)*(-1:2/size(tmp,1):1-1/size(tmp,1)));
        tmp = tmp .* repmat(ph_ramp_t.',[1 nro nrcvrs]);
    end
    
    % Use circular shifts to isolate the desired elements of k
    for r = 1:nrcvrs
        for n=1:Il
            A_t(:,(r-1)*Il+n) =...
                reshape( circshift(tmp(:,:,r),[-I(n),-J(n)]),[],1);
        end
    end
    
    % Compute cost function for all weighting factors
    for wad = 1:szw4
        U_t = A_t*w(:,:,1,wad);
        U_t = reshape(U_t,[npe nro nrcvrs]);
        tester = zeros([npe*2 nro nrcvrs]);
        tester(1:2:end,:,:) = tmp;
        tester(2:2:end,:,:) = U_t;
        
%         if epi_ge == 1
%             % Shift image back
%             ph_ramp_t = exp(-sqrt(-1)*pi*pix(n2)*(-1:2/size(tester,1):1-1/size(tester,1)));
%             tester = tester .* repmat(ph_ramp_t.',[1 nro nrcvrs]);
%         end
        
        tester = fft2_m(tester);
        tester = abs(sqrt(mean(tester.^2,3)));
        %         tester = upsample(tmp(:,:,1),2) + upsample(U_t(:,:,1),2,1);
        %         tester = abs(fftshift(fft2(fftshift(tester))));
        %         ent(wad,n2) = entropy(double(tester/max(max(tester))));

%         ent(wad,n2) = costfcn1(tester) + costfcn2(tester);
%         ent(wad,n2) = costfcn1(tester);
        ent(wad,n2) = cstfcn(tester,low_int_mask);
    end
    
    if all_flag == 1
        % This is for debugging
        tester_all(:,:,n2) = tester;
    end
end

% Find best weighting factors based on entropy 
[~,b_ind_sl,b_shft_sl] = min2D(ent);

end

function out = costfcn1(tester,low_int_mask)
    % Try a cost function that finds the total intensity of the outer
    % regions of the image
    filt =  1-butter2_CB([size(tester,1),1],size(tester,1)*[0.5,10],10);
    filt = repmat(filt,[1 size(tester,2)]);
    tester = tester .* filt / max(max(tester));
    out = sum(sum(tester))/sum(sum(filt));
end

function out = costfcn2(tester,low_int_mask)
    % Try a cost function that finds the total size of a thresholded image
    fit_msk = findsignal(tester,0);
    out = sum(sum(fit_msk))/numel(fit_msk);
end

function out = costfcn3(tester,low_int_mask)
    % Cost function that finds the total intensity of the background
    % voxels, which are identified using the fully sampled reference scan.
    %     fit_msk = logical(round(1-low_int_mask(:,:,1,1,1)));
    out = mean(tester(low_int_mask));
end






