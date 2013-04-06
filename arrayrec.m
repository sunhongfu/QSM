function img = arrayrec(img,fc)
%   Optimal image combination for array coil imaging
%   
%   Usage: img2 = arrayrec(img,f)
%   Author: R Marc Lebel
%       Reference: MRM 47:539-548 (2002)
%   Date: 04/2007
%   
%   Inputs:
%   img: individual array coil complex images np x nv (x ns) x nc.
%   fc:  filter cutoff frequency (def. 0.1)
%   
%   Output:
%   img2: Reconstructed images np x nv (x ns)

%   Check inputs
if nargin < 1
    error('arrayrec: function requires at least one input');
end
if nargin < 2 || isempty(fc)
    fc = 0.1;
end

%   Check image size
np = size(img);
if length(np) < 3 || length(np) > 4
    error('arrayrec: image input must be of size np x nv x (x ns) x nc');
end
if length(np) == 3
    img = reshape(img,[np(1) np(2) 1 np(3)]);
end
[np,nv,ns,nc] = size(img);

%   Check image data type
if isa(img,'double')
    dt = 0;
elseif isa(img,'single')
    dt = 1;
else
    error('arrayrec: image must be of type double or single');
end

%   Create low resolution (filtered) image
%   Convert to k-space
imgl = fftshift(fft2(single(img)));
if dt
    H = butter2(single(ones(np,nv)),fc);
else
    H = butter2(ones(np,nv),fc);
end
imgl = imgl.*repmat(H,[1 1 ns nc]);
clear H
imgl = ifft2(ifftshift(imgl));

%   Create reference and sensitivity images
%   Could improve somewhat using polynomial fitting for b
alpha = sqrt(sum(abs(imgl).^2,4));
b = imgl./repmat(alpha,[1 1 1 nc]);
clear alpha imgl

%   Compute optimal reconstruction
img = 1./sqrt(sum(abs(b).^2,4)) .* sum(img.*conj(b),4);
clear b

%   Retain only real component,and set negative noise to zero
%   May need to apply phase correction term as in paper...
img = real(img);
img(img < 0) = 0;
