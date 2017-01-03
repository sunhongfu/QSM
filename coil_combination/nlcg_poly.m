function [m,RES] = nlcg_poly(m0,params)
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

% iterations
while(k <= params.Itnlim)
   
    % backtracking line-search
    t = t0;
    f0 = objFunc(m,dm,0,params);
    f1 = objFunc(m,dm,t,params);
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
     
    % outputs for debugging purpose
    fprintf('%d , relative changes: %f\n',...
            k, norm(t*dm(:))/norm(m(:)));
    
    if (norm(t*dm(:))/norm(m(:)) <= gradToll); 
        count = count + 1;
    else
        count = 0;
    end 
    if (count == 10)
        break;
    end

    f = f1;
end

return;



function obj = objFunc(m, dm, t, params)
w1 = m+t*dm;
RES = exp(1j*params.P*w1) - params.data;
obj = RES(:)'*RES(:) + params.lambda*w1'*w1;





function grad = wGradient(m,params)

% construct the diagonal matrix

grad = (-1j*params.P.'.*repmat(exp(-1j*params.P*m).',[size(params.P,2),1]))*(exp(1j*params.P*m)-params.data) + (1j*params.P.'.*repmat(exp(1j*params.P*m).',[size(params.P,2),1]))*conj(exp(1j*params.P*m)-params.data) + 2*params.lambda*m;



