function [lfs, res] = echofit(ph, mag, par)
%ECHOFIT Magnitude-weighted least square regression of phase to echo time.
%   [LFS,RES] = echofit(PH,MAG,PAR) fits the phases with TEs, weighted by
%   magnitudes and force intercepts to zeros
%
%   PH:  unwrapped phases from multiple echoes
%   MAG: corresponding magnitudes to phases
%   PAR: parameter sets of the sequence
%   LFS: local field shift after fitting
%   RES: fitting residuals


% check ph and mag have same dimensions
if ~ size(ph)==size(mag)
    error('Input phase and magnitude must be in size');
end

[np,nv,ns,ne] = size(ph);
TE = par.te + (0:ne-1)*par.esp;

ph = permute(ph,[4 1 2 3]);
mag = permute(mag,[4 1 2 3]);


ph = reshape(ph,ne,[]);
mag = reshape(mag,ne,[]);
TE = repmat(TE',[1 np*nv*ns]);

lfs = sum(mag.*ph.*TE,1)./(sum(mag.*TE.*TE)+eps);
lfs = reshape(lfs,[np nv ns]);

lfs_dup = permute(repmat(lfs(:),[1 ne]),[2 1]);
res = -reshape(sum((ph - lfs_dup.*TE).*mag.*(ph - lfs_dup.*TE),1)./sum(mag,1)*ne,[np nv ns]);% % normalize lfs ("/sum*ne")
res(isnan(res)) = 0;
res(isinf(res)) = 0;

