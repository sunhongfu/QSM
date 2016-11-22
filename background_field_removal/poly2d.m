function polyfit = poly2d(lfs,mask,poly_order)


if ~ exist('poly_order','var') || isempty(poly_order)
    poly_order = 2;
end

[np nv ns] = size(lfs);
polyfit = zeros(size(lfs));

% polyfit

px = repmat((1:np)',[nv,1]);
py = repmat((1:nv),[np,1]);
py = py(:);

for i = 1:ns


	lfs_i = lfs(:,:,i);
	mask_i = mask(:,:,i);

	% fit only the non-zero region
	px_nz = px(logical(mask_i(:)));
	py_nz = py(logical(mask_i(:)));


if poly_order == 1
	% first order polyfit
	P = [px, py, ones(length(px),1)]; % polynomials
	P_nz = [px_nz, py_nz, ones(length(px_nz),1)]; % polynomials

elseif poly_order == 2
	% second order polyfit
	P = [px.^2, py.^2, px.*py, px, py, ones(length(px),1)]; % polynomials
	P_nz = [px_nz.^2, py_nz.^2, px_nz.*py_nz, px_nz, py_nz, ones(length(px_nz),1)]; % polynomials

elseif poly_order == 3
	% third order
	P = [px.^3, py.^3, px.*py.^2, py.*px.^2, px.^2, py.^2, px.*py, px, py, ones(length(px),1)]; % polynomials
	P_nz = [px_nz.^3, py_nz.^3, px_nz.*py_nz.^2, py_nz.*px_nz.^2, px_nz.^2, py_nz.^2, px_nz.*py_nz, px_nz, py_nz, ones(length(px_nz),1)]; % polynomials

else
	error('cannot do higher than 3rd order');
end


	I = lfs_i(logical(mask_i)); % measurements of non-zero region
	I = I(:);
%	coeff = P_nz\I; % polynomial coefficients
	coeff = (P_nz'*P_nz)\(P_nz'*I);
	polyfit_i = P*coeff; % residual after polyfit
	polyfit_i = reshape(polyfit_i,[np,nv]);
	polyfit(:,:,i) = polyfit_i;

end


polyfit(isnan(polyfit)) = 0;
polyfit(isinf(polyfit)) = 0;
