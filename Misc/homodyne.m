function [u] = homodyne(u, n, wname, sigma)
%HOMODYNE [u] = homodyne(u, n, [wname, sigma])

% author: from QSM.m github, by Kames UBC

    narginchk(2, 4)

    if nargin < 4, sigma = [0.25, 0.25]; end
    if nargin < 3 || isempty(wname), wname = 'hann'; end

    wname = validatestring(wname, {'hann', 'hamming', 'gaussian'}, 3);

    if length(n) == 1, n = [n, n]; end
    if length(sigma) == 1, sigma = [sigma, sigma]; end

    if n(1) > size(u, 1), n = [size(u, 1), n(2)]; end
    if n(2) > size(u, 2), n = [n(1), size(u, 2)]; end


    w = window_(size(u), n, wname, sigma);

    for t = 1:size(u, 4)
        for k = 1:size(u, 3)
            slice = u(:,:,k,t);
            u(:,:,k,t) = slice ./ ifft2(w .* fft2(slice));
        end
    end

end


function [w] = window_(sz, n, wname, sigma)

    n1 = n(1);
    n2 = n(2);
    n12 = floor(n1 / 2);
    n22 = floor(n2 / 2);

    x = [linspace(0, n12/n1, n12+1), linspace(-n12/n1, -1/n1, n12)].';
    y = [linspace(0, n22/n2, n22+1), linspace(-n22/n2, -1/n2, n22)].';

    switch wname
        case 'hann'
            x = hann_(x);
            y = hann_(y);
        case 'hamming'
            x = hamming_(x);
            y = hamming_(y);
        case 'gaussian'
            x = gaussian_(x, sigma(1));
            y = gaussian_(y, sigma(2));
    end

    wx = zeros(sz(1), 1);
    wy = zeros(sz(2), 1);

    wx(1:n12+1) = x(1:n12+1);
    wx(end-n12+1:end) = x(n12+2:end);

    wy(1:n22+1) = y(1:n22+1);
    wy(end-n22+1:end) = y(n22+2:end);

    w = wx * wy.';

end


function [x] = hann_(x)
    x = 0.5 .* (1 + cos(2*pi .* x));
end


function [x] = hamming_(x)
    x = 0.54 + 0.46.*cos(2*pi .* x);
end


function [x] = gaussian_(x, sigma)
    x = exp(-0.5 .* (x./sigma).^2);
end
