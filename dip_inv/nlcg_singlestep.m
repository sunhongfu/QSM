function [m,Res_term,TV_term] = nlcg_singlestep(m0,params)
% Phi(m) = ||W(Fu*m - y)||^2 + lamda1*|TV*m|_1
% m: susceptibility
% W: weighting matrix derived from magnitude intensities
% Fu: F_{-1}*D*F forward calculates the field from susceptibility
% y: measured field to be fitted (inversion)
% lambda1: TV regularization parameter
% TV: total variation operation
% ||...||^2: L2 norm
% |...|_1: L1 norm
% note the TV term can also be L2 norm if set p=2,
% then the term would be changed to ||TV*m||^2



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

f = 0;
Res_term0 = 0;
count =0;

% iterations
while(k <= params.Itnlim)
   
    % backtracking line-search
    t = t0;
    f0 = objFunc(m,dm,0,params);
    [f1, Res_term, TV_term] = objFunc(m,dm,t,params);
    lsiter = 0;
    while (f1 > f0 + alpha*t*(g0(:)'*dm(:))) && (lsiter<maxlsiter)
        t = t* beta;
        [f1, Res_term, TV_term] = objFunc(m,dm,t,params);
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
     
    % outputs for debugging purpose
    % fprintf('%d , relative residual: %f\n',...
    %         k, abs(Res_term-Res_term0)/Res_term);
    
    % if (abs(Res_term-Res_term0)/Res_term <= gradToll); 
    %     count = count + 1;
    % else
    %     count = 0;
    % end 
    
    % fprintf('%d , relative changes: %f\n',...
    %         k, norm(t*dm(:))/norm(m(:)));
    fprintf('.');
    
    if (norm(t*dm(:))/norm(m(:)) <= gradToll); 
        count = count + 1;
    else
        count = 0;
    end 
    if (count == 10)
        break;
    end

    f = f1;
    Res_term0 = Res_term;
end

return;



function [obj,Res_term,TV_term,Tik_term] = objFunc(m, dm, t, params)
p = params.pNorm;
w1 = m+t*dm;

w2 = params.TV*(w1.*params.TV_mask);
TV = (w2.*conj(w2)+params.l1Smooth).^(p/2);
TV_term = sum(params.TV_reg(:).*TV(:));

Res_term = params.FT*(w1.*params.sus_mask) - params.data;
Res_term = (params.Res_wt(:).*Res_term(:))'*(params.Res_wt(:).*Res_term(:));

Tik_term = params.Tik_reg*((params.Tik_mask(:).*w1(:))'*((params.Tik_mask(:).*w1(:))));

obj = Res_term + Tik_term + TV_term;



function grad = wGradient(m,params)
p = params.pNorm;
w1 = params.TV*(m.*params.TV_mask);

grad_TV = params.TV_mask.*(params.TV'*(p*w1.*(w1.*conj(w1)+params.l1Smooth).^(p/2-1)));

grad_Res = params.sus_mask.*(params.FT'*((params.Res_wt.^2).*((params.FT*(params.sus_mask.*m))-params.data)));

grad_Tik = params.Tik_mask.*m;

grad = 2*grad_Res + params.TV_reg*grad_TV + 2*params.Tik_reg.*grad_Tik;

