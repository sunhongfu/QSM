function [lfs, res] = echofit(ph, mag, te)
%ECHOFIT Magnitude-weighted least square regression of phase to echo time.
%   [LFS,RES] = echofit(PH,MAG,PAR) fits the phases with TEs, weighted by
%   magnitudes and force intercepts to zeros
%
%   PH:  unwrapped phases from multiple echoes
%   MAG: corresponding magnitudes to phases
%   TE : echo times in a row (not column), e.g. [3 7 11 15 19] (ms)
%   LFS: local field shift after fitting
%   RES: fitting residuals


% check ph and mag have same dimensions
if ~ size(ph)==size(mag)
    error('Input phase and magnitude must be in size');
end

[np,nv,ns,ne] = size(ph);

ph = permute(ph,[4 1 2 3]);
mag = permute(mag,[4 1 2 3]);


ph = reshape(ph,ne,[]);
mag = reshape(mag,ne,[]);
te = repmat(te',[1 np*nv*ns]);

lfs = sum(mag.*ph.*te,1)./(sum(mag.*te.*te)+eps);
lfs = reshape(lfs,[np nv ns]);

lfs_dup = permute(repmat(lfs(:),[1 ne]),[2 1]);
res = reshape(sum((ph - lfs_dup.*te).*mag.*(ph - lfs_dup.*te),1)./sum(mag,1)*ne,[np nv ns]);% % normalize lfs ("/sum*ne")
res(isnan(res)) = 0;
res(isinf(res)) = 0;

