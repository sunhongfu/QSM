function Mh = butter2(M,kc)
%   Applies a low pass Butterworth filter to the frequency data
%
%   Author: Marc Lebel 07/06
%   Usage:
%   Mh = butter2(M,kc)
%   To output the filter, input a unit matrix
%   
%   input:
%   M is a 1, 2, or 3 dimensional matrix
%   kc is the cuttoff frequency. Default half
%   
%   OUTPUT:
%   Mh is the filtered data of the same size

%   Check base input argument
if nargin < 1
    error('butter2 requires at least one input')
end

%   Check k space size
N = size(M);
if length(N) < 2 || length(N) > 3
    error('butter2: wrong data size');
end

%   Assign default cutoff frequency
if nargin < 2
    kc = 0.7;
end

%   Create 1D filter
if N(1) == 1 || N(2) == 1
    k1 = -1:2/(length(M)-1):1;
    B = (1./(1+((k1./kc).^2).^4));
elseif ndims(M) == 2
    %   Create 2D filter
    [k1,k2] = meshgrid(-1:2/(N(2)-1):1,-1:2/(N(1)-1):1);
    L = k1.^2 + k2.^2;
    B = 1./(1+(L/kc.^2).^4);
    clear k1 k2 L
elseif ndims(M) == 3
    [k1,k2,k3] = ndgrid(-1:2/(N(1)-1):1,-1:2/(N(2)-1):1,-1:2/(N(3)-1):1);
    L = k1.^2 + k2.^2 + k3.^2;
    B = 1./(1+(L/kc.^2).^4);
end

%   Apply filter
Mh = M.*B;
