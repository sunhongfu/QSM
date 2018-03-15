function [m, Res_term, TV_term, Tik_term] = nlcg_tik_gpu(m0, params)
%
% m: susceptibility chi of whole FOV or only local brain tissue
% Res_term: norm of residual/fidelity term (L-curve purpose)
% TV_term: norm of TV term
% Tik_term: norm of Tik_term
%
% argmin ||Res_wt * (F_{-1} * D * F * sus_mask * chi - tfs)|| + TV_reg * TV|TV_mask * chi| + Tik_reg * ||Tik_mask * chi|| + TV_reg2 * TV|air_mask * chi|
%
% tfs:      total field shift; can be local field shift if use this for local field inversion
% Res_wt:   weighting matrix for the residual/fidelity term, usually brain mask
% sus_mask: input 1 for total field inversion; input brain mask for local field inversion
% TV_mask:  usually brain mask for TV regularization of the local tissue susceptiblity distribution
% Tik_mask: usually brain mask for Tikhonov regularization of the local tissue susceptiblity distribution
% TV_reg:   Regularization parameter for TV term
% Tik_reg:  Regularization parameter for Tikhonov term


m = gpuArray(m0);

% line search parameters
maxlsiter = gpuArray(params.lineSearchItnlim);
gradToll  = gpuArray(params.gradToll);
alpha     = gpuArray(params.lineSearchAlpha); 
beta      = gpuArray(params.lineSearchBeta);
t0        = gpuArray(params.lineSearchT0);
k         = gpuArray(0);
Itnlim    = gpuArray(params.Itnlim);


P         = gpuArray(params.P);
TV_mask   = gpuArray(params.TV_mask);
sus_mask  = gpuArray(params.sus_mask);
Res_wt    = gpuArray(params.Res_wt);
Tik_mask  = gpuArray(params.Tik_mask);
Tik_reg   = gpuArray(params.Tik_reg);
TV_reg2   = gpuArray(TV_reg2);
air_mask  = gpuArray(params.air_mask);
TV_reg    = gpuArray(params.TV_reg);
pNorm     = gpuArray(params.pNorm);
data      = gpuArray(params.data);
l1Smooth  = gpuArray(params.l1Smooth);
D         = gpuArray(params.D);



% compute (-gradient): search direction
g0 = wGradient(m,pNorm, P, TV_mask, air_mask, l1Smooth, sus_mask, D, data, Tik_mask, TV_reg, TV_reg2);
dm = -g0;

% f = 0;
% Res_term0 = 0;
count = gpuArray(0);

% iterations
while(k <= Itnlim)
   
    % backtracking line-search
    t = t0;
    f0 = objFunc(m,dm,0,pNorm, P, TV_mask, l1Smooth, air_mask, sus_mask, D, data, Tik_mask, TV_reg, TV_reg2);
    [f1, Res_term, TV_term, Tik_term, TV_term2] = objFunc(m, dm, t, pNorm, P, TV_mask, l1Smooth, air_mask, sus_mask, D, data, Tik_mask, TV_reg, TV_reg2);
    lsiter = gpuArray(0);
    while (f1 > f0 + alpha*t*(g0(:)'*dm(:))) && (lsiter<maxlsiter)
        t = t* beta;
        f1 = objFunc(m,dm,t,pNorm, P, TV_mask, l1Smooth, air_mask, sus_mask, D, data, Tik_mask, TV_reg, TV_reg2);
        lsiter = lsiter + 1;
    end
    % control the number of line searches by adapting the initial step search
    if lsiter > 2
		t0 = t0 * beta;
    end 
    if lsiter < 1
		t0 = t0 / beta;
    end
    
    % updates
    m = m + t*dm; dm0 = dm;
    g1 = wGradient(m,pNorm, P, TV_mask, air_mask, l1Smooth, sus_mask, D, data, Tik_mask, TV_reg, TV_reg2);
    bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
    g0 = g1;
    dm = -g1 + bk*dm;
    k = k + 1;
     

    fprintf('.');
    
    if (norm(t*dm(:))/norm(m(:)) <= gradToll); 
        count = count + 1;
    else
        count = 0;
    end 
    if (count == 10) % under gradToll for continueous 10 times
        break;
    end


end

return;



function [obj,Res_term,TV_term,Tik_term,TV_term2] = objFunc(m, dm, t, pNorm, P, TV_mask, l1Smooth, air_mask, sus_mask, D, data, Tik_mask, TV_reg, TV_reg2)
p = pNorm;
w1 = m+t*dm;

w2 = Diff(P.*w1.*TV_mask);
TV = (w2.*conj(w2)+l1Smooth).^(p/2);
TV_term = sum(TV(:));

% add air TV
w3 = Diff(P.*w1.*air_mask);
TV2 = (w3.*conj(w3)+l1Smooth).^(p/2);
TV_term2 = sum(TV2(:));

Res_term = ifftn(fftn((P.*w1.*sus_mask).*D)) - data;
Res_term = (Res_wt(:).*Res_term(:))'*(Res_wt(:).*Res_term(:));

Tik_term = (P(:).*Tik_mask(:).*w1(:))'*(P(:).*Tik_mask(:).*w1(:));

% obj = Res_term + Tik_reg*Tik_term + TV_reg*TV_term;
obj = Res_term + Tik_reg*Tik_term + TV_reg*TV_term + TV_reg2*TV_term2;



function grad = wGradient(m, pNorm, P, TV_mask, air_mask, l1Smooth, sus_mask, D, data, Tik_mask, TV_reg, TV_reg2)
p = pNorm;
w1 = Diff(P.*m.*TV_mask);
w2 = Diff(P.*m.*air_mask);

grad_TV = P.*TV_mask.*(invD(p*w1.*(w1.*conj(w1)+l1Smooth).^(p/2-1)));
grad_TV2 = P.*air_mask.*(invD(p*w2.*(w2.*conj(w2)+l1Smooth).^(p/2-1)));

grad_Res = P.*sus_mask.*ifftn(D.*fftn((Res_wt.^2).*(ifftn(fftn((P.*sus_mask.*m).*D))-data)));

grad_Tik = P.^2.*Tik_mask.^2.*m;

grad = 2*grad_Res + TV_reg*grad_TV + 2*Tik_reg.*grad_Tik + TV_reg2*grad_TV2;





function res = invD(y)
res = adjDx(y(:,:,:,1)) + adjDy(y(:,:,:,2)) + adjDz(y(:,:,:,3));

function res = Diff(x)
Dx = x([2:end,end],:,:) - x;
Dy = x(:,[2:end,end],:) - x;
Dz = x(:,:,[2:end,end]) - x;
res = cat(4,Dx,Dy,Dz);

function res = adjDx(x)
res = x([1,1:end-1],:,:) - x;
res(1,:,:) = -x(1,:,:);
res(end,:,:) = x(end-1,:,:);

function res = adjDy(x)
res = x(:,[1,1:end-1],:) - x;
res(:,1,:) = -x(:,1,:);
res(:,end,:) = x(:,end-1,:);

function res = adjDz(x)
res = x(:,:,[1,1:end-1]) - x;
res(:,:,1) = -x(:,:,1);
res(:,:,end) = x(:,:,end-1);
