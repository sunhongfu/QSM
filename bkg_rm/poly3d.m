function polyfit = poly3d(lfs,mask)

[np nv nv2] = size(lfs);

% polyfit
px = repmat((1:np)',[nv*nv2,1]);
py = repmat((1:nv),[np,1]);
py = repmat(py(:),[nv2,1]);
pz = repmat((1:nv2),[np*nv,1]);
pz = pz(:);

% fit only the non-zero region
px = px(logical(mask(:)));
py = py(logical(mask(:)));
pz = pz(logical(mask(:)));

% first order polyfit
% P = [px, py, pz, ones(length(px),1)]; % polynomials

% second order
P = [px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials

I = lfs(logical(mask)); % measurements of non-zero region
I = I(:);
coeff = P\I; % polynomial coefficients
residual = I - P*coeff; % residual after polyfit

% name the phase result after polyfit as tfs (total field shift)
polyfit = zeros(np*nv*nv2,1);
polyfit(logical(mask(:))) = residual;
polyfit = reshape(polyfit,[np,nv,nv2]);
