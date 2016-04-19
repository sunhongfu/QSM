function polyfit = poly3d(lfs,mask,poly_order)

if ~ exist('poly_order','var') || isempty(poly_order)
    poly_order = 2;
end

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

if poly_order == 1
% first order polyfit
P = [px, py, pz, ones(length(px),1)]; % polynomials

elseif poly_order == 2	
% second order
P = [px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials

elseif poly_order == 3
% third order
P = [px.^3, py.^3, pz.^3, px.*py.^2, px.*pz.^2, px.*py.*pz, py.*px.^2, py.*pz.^2, pz.*px.^2, pz.*py.^2 ...
	px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials

else
	error('cannot do higher than 3rd order');
end


I = lfs(logical(mask)); % measurements of non-zero region
I = I(:);
% coeff = P\I; % polynomial coefficients
coeff = (P'*P)\(P'*I);
residual = I - P*coeff; % residual after polyfit

% name the phase result after polyfit as tfs (total field shift)
polyfit = zeros(np*nv*nv2,1);
polyfit(logical(mask(:))) = residual;
polyfit = reshape(polyfit,[np,nv,nv2]);


polyfit(isnan(polyfit)) = 0;
polyfit(isinf(polyfit)) = 0;
