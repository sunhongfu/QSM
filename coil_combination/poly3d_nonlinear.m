% nonlinear polynomial fit using nonlinear CG 
function polyfit_nonlinear = poly3d_nonlinear(offsets,mask,poly_order)

if ~ exist('poly_order','var') || isempty(poly_order)
    poly_order = 1;
end

[np nv nv2] = size(offsets);

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
end

I = offsets(logical(mask)); % measurements of non-zero region
I = I(:);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start non-linear CG method
params.Itnlim = 200; % interations numbers (adjust accordingly!)
params.gradToll = 1e-4; % step size tolerance stopping criterea
params.lineSearchItnlim = 200;
params.lineSearchAlpha = 0.01;
params.lineSearchBeta = 0.6;
params.lineSearchT0 = 1 ; % step size to start with
params.data = I;
params.P = P;
% params.wt = 1; % weighting matrix
params.lambda = 1e6;



coeff = nlcg_poly(zeros(4,1), params);

% name the phase result after polyfit as tfs (total field shift)
polyfit_nonlinear = zeros(np*nv*nv2,1);
polyfit_nonlinear(logical(mask(:))) = exp(1j*P*coeff);
polyfit_nonlinear = reshape(polyfit_nonlinear,[np,nv,nv2]);

nii = make_nii(angle(polyfit_nonlinear));
save_nii(nii,'offsets_poly.nii')
