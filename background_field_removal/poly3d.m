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
px_nz = px(logical(mask(:)));
py_nz = py(logical(mask(:)));
pz_nz = pz(logical(mask(:)));

if poly_order == 1
	% first order polyfit
	P = [px, py, pz, ones(length(px),1)];
	P_nz = [px_nz, py_nz, pz_nz, ones(length(px_nz),1)]; % polynomials

elseif poly_order == 2	
	% second order
	P = [px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials
	P_nz = [px_nz.^2, py_nz.^2, pz_nz.^2, px_nz.*py_nz, px_nz.*pz_nz, py_nz.*pz_nz, px_nz, py_nz, pz_nz, ones(length(px_nz),1)]; % polynomials

elseif poly_order == 3
	% third order
	P = [px.^3, py.^3, pz.^3, px.*py.^2, px.*pz.^2, px.*py.*pz, py.*px.^2, py.*pz.^2, pz.*px.^2, pz.*py.^2 ...
		px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials
	P_nz = [px_nz.^3, py_nz.^3, pz_nz.^3, px_nz.*py_nz.^2, px_nz.*pz_nz.^2, px_nz.*py_nz.*pz_nz, py_nz.*px_nz.^2, py_nz.*pz_nz.^2, pz_nz.*px_nz.^2, pz_nz.*py_nz.^2 ...
		px_nz.^2, py_nz.^2, pz_nz.^2, px_nz.*py_nz, px_nz.*pz_nz, py_nz.*pz_nz, px_nz, py_nz, pz_nz, ones(length(px_nz),1)]; % polynomials


elseif poly_order == 4
	% third order
	P = [px.^4, py.^4, pz.^4, px.^3.*py, px.^3.*pz, py.^3.*px, py.^3.*pz, pz.^3.*px, pz.^3.*py, px.^2.*py.^2, px.^2.*pz.^2, px.^2.*py.*pz, py.^2.*pz.^2, py.^2.*px.*pz, pz.^2.*px.*py ...
		px.^3, py.^3, pz.^3, px.*py.^2, px.*pz.^2, px.*py.*pz, py.*px.^2, py.*pz.^2, pz.*px.^2, pz.*py.^2 ...
		px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials
	P_nz = [px_nz.^4, py_nz.^4, pz_nz.^4, px_nz.^3.*py_nz, px_nz.^3.*pz_nz, py_nz.^3.*px_nz, py_nz.^3.*pz_nz, pz_nz.^3.*px_nz, pz_nz.^3.*py_nz, px_nz.^2.*py_nz.^2, px_nz.^2.*pz_nz.^2, px_nz.^2.*py_nz.*pz_nz, py_nz.^2.*pz_nz.^2, py_nz.^2.*px_nz.*pz_nz, pz_nz.^2.*px_nz.*py_nz ...
		px_nz.^3, py_nz.^3, pz_nz.^3, px_nz.*py_nz.^2, px_nz.*pz_nz.^2, px_nz.*py_nz.*pz_nz, py_nz.*px_nz.^2, py_nz.*pz_nz.^2, pz_nz.*px_nz.^2, pz_nz.*py_nz.^2 ...
		px_nz.^2, py_nz.^2, pz_nz.^2, px_nz.*py_nz, px_nz.*pz_nz, py_nz.*pz_nz, px_nz, py_nz, pz_nz, ones(length(px_nz),1)]; % polynomials

else
	error('cannot do higher than 3rd order');
end


I = lfs(logical(mask)); % measurements of non-zero region
I = I(:);
% coeff = P_nz\I; % polynomial coefficients
coeff = (P_nz'*P_nz)\(P_nz'*I);
polyfit = P*coeff;
polyfit = reshape(polyfit,[np,nv,nv2]);


polyfit(isnan(polyfit)) = 0;
polyfit(isinf(polyfit)) = 0;
