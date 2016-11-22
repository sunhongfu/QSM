function w = grappa_findw(k,opt)
%   Determines GRAPPA weightings for an arbitrary kernal. This should
%   support any undersampling rate (which is reflected in the kernal);
%   however, it has only been tested for rate 2 undersampling. 
%   
%   Usage: w = grappa_findw(k,opt)
%   Author: Corey Baron
%   Date: 10/2010
%   
%   Input:
%   k: k-space data for fully sampled reference data of size
%       Npe x Nro x NS x (# calibration scans) x RCVRS
%       This data should already be Nyquist ghost corrected (for EPI) and gridded.
%   opt:    structure with recon options
%       opt.kern: kernal used for GRAPPA weights. A value of -1 corresponds to the
%           center of the kernal, 0 for unused points, and 1 for used points. The
%           kernal must include lines that will be missing in the undersampled data
%           (these should be rows of zeros in the kernal).
%       opt.kernreg: regulization factor. 10^-2 by default.
%       opt.PI_multw: When multiple calibration scans are provided, 
%           all calibration scans can be used together in matrix inversion 
%           (opt.PI_multw = 0) (recommended), or seperate weighting factors 
%           can be determined for each reference scan (opt.PI_multw = 1).
%       opt.numwork: number of workers to use for parallel processing.
%   
%   Output:
%   w: GRAPPA weights of size 
%       (RCVRS * # non-zero kernal entries) x RCVRS x NS)

% load refdata (for testing)
% k = squeeze(k_ref(:,:,:,end,:));

% Choose whether to exclude edge voxels in determining the the wieghts
% (since I use circular shifts, they should be excluded)
exclude_edgevoxels = 1;

% Move the slices and # of calibration scans to the end
k = permute(k,[1 2 5 4 3]);
%k = permute(k,[1 2 4 3]);
npnts = size(k,1)*size(k,2);
nrcvrs = size(k,3);
nslices = size(k,5);
nwad = size(k,4);

% Check options
if ~isfield(opt,'kernreg') || isempty(opt.kernreg)
    opt.kernreg = 10^-2;
end
if ~isfield(opt,'PI_multw') || isempty(opt.PI_multw)
    opt.PI_multw = 0;
end
if ~isfield(opt,'numwork') || isempty(opt.numwork)
    opt.numwork = 1;
end

% NB: in the kernal, 0 corresponds to values not used, 1 to values used, -1
% to the center.
[Ci, Cj] = find(opt.kern == -1);
if length(Ci) > 1
    error('Grappa kernal inadequately defined.')
end

% Find the non-zero entries in the kernal
[I,J] = find(opt.kern>0);
I = I - Ci;
J = J - Cj;
mI = max(abs(I));
mJ = max(abs(J));
npnts_A = npnts;
if exclude_edgevoxels == 1
    extra_ex = 0;
    mI = mI+extra_ex;
    mJ = mJ+extra_ex;
    npnts_A = (size(k,1)-2*mI)*(size(k,2)-2*mJ);
end

% Perform the calculations (without regularization)
% U = reshape(k,[npnts nrcvrs nslices nwad]);
% A = zeros(npnts,length(I)*nrcvrs,nslices, nwad);
% w = zeros(length(I)*nrcvrs,nrcvrs,nslices, nwad);
% for s = 1:nslices*nwad
%     % Use circular shifts to isolate the desired elements of k
%     for r = 1:nrcvrs
%         for n=1:length(I)
%             A(:,(r-1)*length(I)+n,s) =...
%                 reshape( circshift(k(:,:,r,s),[-I(n),-J(n)]),[],1);
%         end
%     end
%     
%     % Calculate the weights
%     w(:,:,s) = A(:,:,s)\U(:,:,s);
% end

% Convert images to magnitudes
    % NB: just thought I'd try it: it had quite disastrous effects
% sz = size(k);
% for n=1:prod(sz(3:end))
%     k(:,:,n) = fftshift(fft2(fftshift(squeeze(k(:,:,n)))));
%     k(:,:,n) = ifftshift(ifft2(ifftshift(squeeze(abs(k(:,:,n))))));
% end

if opt.PI_multw == 1
    nwad2 = nwad;
    nwad = 1;
    w = zeros(length(I)*nrcvrs,nrcvrs,nslices,nwad2)+1i*ones(length(I)*nrcvrs,nrcvrs,nslices,nwad2);
else
    nwad2 = 1;
    w = zeros(length(I)*nrcvrs,nrcvrs,nslices)+1i*ones(length(I)*nrcvrs,nrcvrs,nslices);
end
for m = 1:nwad2
    % Compute using regularization
    if opt.PI_multw == 1
        if exclude_edgevoxels == 1
            U = permute(k(1+mI:end-mI,1+mJ:end-mJ,:,:,:), [1 2 3 5 4]);
        else
            U = permute(k, [1 2 3 5 4]);
        end
        U = reshape(U,[npnts_A nrcvrs nslices nwad2]);
    else
        if exclude_edgevoxels == 1
            U = permute(k(1+mI:end-mI,1+mJ:end-mJ,:,:,:), [1 2 4 3 5]);
        else
            U = permute(k, [1 2 4 3 5]);
        end
        U = reshape(U,[npnts_A*nwad nrcvrs nslices]);
    end
    if opt.numwork > 0
        % Use parallel processing if available
        U2 = U(:,:,:,m);
        if opt.PI_multw == 1
            k2 = k(:,:,:,m,:);
        else
            k2 = k;
        end
        PI_multw = opt.PI_multw;
        parfor (s = 1:nslices,opt.numwork)
            %     for s = 1:nslices
            A = zeros(npnts_A*nwad,length(I)*nrcvrs)+1i*ones(npnts_A*nwad,length(I)*nrcvrs);
            
            % Use circular shifts to isolate the desired elements of k
            k_a = k2(:,:,:,:,s);
            for a=1:nwad
                for r = 1:nrcvrs
                    for n=1:length(I)
                        if exclude_edgevoxels == 1
                            shifted_data = circshift(k_a(:,:,r,a),[-I(n),-J(n)]);
                            shifted_data = shifted_data(1+mI:end-mI,1+mJ:end-mJ,:,:);
                            A( (1:npnts_A)+(a-1)*npnts_A,(r-1)*length(I)+n) =...
                                reshape( shifted_data,[],1);
                        else
                            A( (1:npnts)+(a-1)*npnts,(r-1)*length(I)+n) =...
                                reshape( circshift(k_a(:,:,r,a),[-I(n),-J(n)]),[],1);
                        end
                    end
                end
            end
            
            % Prepare for inversion
            At = A'*A;
            lambda = norm(At,'fro')/size(At,1) * opt.kernreg;
            
            % Calculate the weights
            w_a = zeros(length(I)*nrcvrs,nrcvrs)+1i*ones(length(I)*nrcvrs,nrcvrs);
            U_a = U2(:,:,s);
            for rcvr = 1:nrcvrs
                w_a(:,rcvr) = pinv(At + lambda*eye(size(At))) * A' * U_a(:,rcvr);
            end
            w(:,:,s,m) = w_a;
        end
        clear U2 k2
    else
        A = zeros(npnts_A*nwad,length(I)*nrcvrs)+1i*ones(npnts_A*nwad,length(I)*nrcvrs);
        for s = 1:nslices
            % Use circular shifts to isolate the desired elements of k
            if opt.PI_multw == 1
                k_a = k(:,:,:,m,s);
            else
                k_a = k(:,:,:,:,s);
            end
            for a=1:nwad
                for r = 1:nrcvrs
                    for n=1:length(I)
                        if exclude_edgevoxels == 1
                            shifted_data = circshift(k_a(:,:,r,a),[-I(n),-J(n)]);
                            shifted_data = shifted_data(1+mI:end-mI,1+mJ:end-mJ,:,:);
                            A( (1:npnts_A)+(a-1)*npnts_A,(r-1)*length(I)+n) =...
                                reshape( shifted_data,[],1);
                        else
                            A( (1:npnts)+(a-1)*npnts,(r-1)*length(I)+n) =...
                                reshape( circshift(k_a(:,:,r,a),[-I(n),-J(n)]),[],1);
                        end
                    end
                end
            end
            
            % Prepare for inversion
            At = A'*A;
            lambda = norm(At,'fro')/size(At,1) * opt.kernreg;
            
            % Calculate the weights
            w_a = zeros(length(I)*nrcvrs,nrcvrs)+1i*ones(length(I)*nrcvrs,nrcvrs);
            U_a = U(:,:,s,m);
            for rcvr = 1:nrcvrs
                w_a(:,rcvr) = pinv(At + lambda*eye(size(At))) * A' * U_a(:,rcvr);
            end
            w(:,:,s,m) = w_a;
        end
    end
    
end

end
