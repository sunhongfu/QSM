function [m, Res_term, TV_term, Tik_term] = nlcg_tik(m0, params)
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


m = m0;

% line search parameters
maxlsiter = params.lineSearchItnlim;
gradToll  = params.gradToll;
alpha     = params.lineSearchAlpha; 
beta      = params.lineSearchBeta;
t0        = params.lineSearchT0;
k         = 0;

% compute (-gradient): search direction
g0 = wGradient(m,params);
dm = -g0;

% f = 0;
% Res_term0 = 0;
count =0;

% iterations
while(k <= params.Itnlim)
   
    % backtracking line-search
    t = t0;
    f0 = objFunc(m,dm,0,params);
    [f1, Res_term, TV_term, Tik_term, TV_term2] = objFunc(m,dm,t,params);
    lsiter = 0;
    while (f1 > f0 + alpha*t*(g0(:)'*dm(:))) && (lsiter<maxlsiter)
        t = t* beta;
        f1 = objFunc(m,dm,t,params);
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
    g1 = wGradient(m,params);
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



function [obj,Res_term,TV_term,Tik_term,TV_term2] = objFunc(m, dm, t, params)
p = params.pNorm;
w1 = m+t*dm;

w2 = params.TV*(params.P.*w1.*params.TV_mask);
TV = (w2.*conj(w2)+params.l1Smooth).^(p/2);
TV_term = sum(TV(:));

% add air TV
w3 = params.TV*(params.P.*w1.*params.air_mask);
TV2 = (w3.*conj(w3)+params.l1Smooth).^(p/2);
TV_term2 = sum(TV2(:));

Res_term = params.FT*(params.P.*w1.*params.sus_mask) - params.data;
Res_term = (params.Res_wt(:).*Res_term(:))'*(params.Res_wt(:).*Res_term(:));

Tik_term = (params.P(:).*params.Tik_mask(:).*w1(:))'*(params.P(:).*params.Tik_mask(:).*w1(:));

% obj = Res_term + params.Tik_reg*Tik_term + params.TV_reg*TV_term;
obj = Res_term + params.Tik_reg*Tik_term + params.TV_reg*TV_term + params.TV_reg2*TV_term2;



function grad = wGradient(m,params)
p = params.pNorm;
w1 = params.TV*(params.P.*m.*params.TV_mask);
w2 = params.TV*(params.P.*m.*params.air_mask);

grad_TV = params.P.*params.TV_mask.*(params.TV'*(p*w1.*(w1.*conj(w1)+params.l1Smooth).^(p/2-1)));
grad_TV2 = params.P.*params.air_mask.*(params.TV'*(p*w2.*(w2.*conj(w2)+params.l1Smooth).^(p/2-1)));

grad_Res = params.P.*params.sus_mask.*(params.FT'*((params.Res_wt.^2).*((params.FT*(params.P.*params.sus_mask.*m))-params.data)));

grad_Tik = params.P.^2.*params.Tik_mask.^2.*m;

grad = 2*grad_Res + params.TV_reg*grad_TV + 2*params.Tik_reg.*grad_Tik + params.TV_reg2*grad_TV2;

