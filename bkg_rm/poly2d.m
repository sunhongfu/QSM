function polyfit = poly2d(lfs,mask)

[np nv ns] = size(lfs);
polyfit = zeros(size(lfs));


for i = 1:ns

	% polyfit
	px = repmat((1:np)',[nv,1]);
	py = repmat((1:nv),[np,1]);
	py = py(:);


	lfs_i = lfs(:,:,i);
	mask_i = mask(:,:,i);

	% fit only the non-zero region
	px = px(logical(mask_i(:)));
	py = py(logical(mask_i(:)));


	% second order polyfit
	% P = [px, py, ones(length(px),1)]; % polynomials
	P = [px.^2, py.^2, px.*py, px, py, ones(length(px),1)]; % polynomials



	I = lfs_i(logical(mask_i)); % measurements of non-zero region
	I = I(:);
	coeff = P\I; % polynomial coefficients
	residual = I - P*coeff; % residual after polyfit

	polyfit_i = zeros(np*nv,1);
	polyfit_i(logical(mask_i(:))) = residual;
	polyfit_i = reshape(polyfit_i,[np,nv]);

	polyfit(:,:,i) = polyfit_i;

	polyfit(isnan(polyfit)) = 0;
	polyfit(isinf(polyfit)) = 0;
end
