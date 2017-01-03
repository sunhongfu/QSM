function [img,par] = r2s47_recon(PATH_IN)
%R2S47_RECON Reconstruct fid file.
%   [IMG,PAR] = r2s47_recon(PATH_IN)
%
%   PATH_IN : directory contains .fid file
%   IMG     : raw complex images
%   PAR     : parameter sets of the sequence


%   Read parameters
par = readprocpar(PATH_IN);

%   Define sequence variables
np     = par.np/2;
nv     = par.nv;
nv2    = par.nv2;
ne     = par.ne;
rcvrs  = par.nrcvrs;
ad     = par.arraydim/(rcvrs*nv2);
lpe    = par.lpe;
lpe2   = par.lpe2;
ppe    = par.ppe;
ppe2   = par.ppe2;

%   Update par
par.pss = -lpe2/2-lpe2/nv2/2:lpe2/nv2:lpe2/2;
par.pss = par.pss + ppe2;

%   Read data and DC correct
k = readfid(PATH_IN,par,1,0,0);

%   Reshape and permute k
rvec = [np ne nv rcvrs nv2 ad];
pvec = [1 3 5 2 4 6];
k = permute(reshape(k,rvec),pvec);

%   Remove DC
for i = 1:ne*ad*rcvrs
    kt = k(:,:,:,i);
    dc = complex(median(real(kt(:))),median(imag(kt(:))));
    k(:,:,:,i) = k(:,:,:,i) - dc;
end

%   Apply ppe phase ramp
if abs(ppe) > 0.001 
    pix = 0.5 * ppe./(lpe./nv);
    ph_ramp = exp(sqrt(-1)*2*pi*pix*(-1:2/nv:1-1/nv));
    ph_ramp = single(ph_ramp);
    k = k .* repmat(ph_ramp,[np 1 nv2 ne rcvrs ad]);
    clear ph_ramp pix
end

%   Apply ppe2 phase ramp
if abs(ppe2) > 0.001
    pix = 0.5 * ppe2./(lpe2./nv2);
    ph_ramp = exp(sqrt(-1)*2*pi*pix*(-1:2/nv2:1-1/nv2));
    ph_ramp = reshape(ph_ramp,[1 1 nv2]);
    ph_ramp = single(ph_ramp);
    k = k .* repmat(ph_ramp,[np nv 1 ne rcvrs ad]);
    clear ph_ramp pix
end

%   Flip data for proper orientation
k = flipdim(flipdim(flipdim(k,1),2),3);

%   Convert to image space
img = zeros(size(k),'single');
for i=1:ne*rcvrs*ad
    %   Take 3rd dimension transform
    img(:,:,:,i) = fftshift(fft(fftshift(k(:,:,:,i),3),[],3),3);

    %   Further remove baseline for DC suppression
    for j = 1:nv2
        itmp = img(:,:,j,i);
        ind = abs(itmp(:)) < 4 * mean(abs(itmp(:)));
        img(:,:,j,i) = img(:,:,j,i) - mean(itmp(ind));
    end
    
    %   Take inplane transform
    img(:,:,:,i) = fftshift(fftshift(fft(fft(fftshift(fftshift( ...
        img(:,:,:,i),1),2),[],1),[],2),1),2);
end
