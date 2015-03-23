function [m,RES,TVterm] = nlcg(m0,params)
% Phi(m) = ||Fu*m - y||^2 + lamda1*|TV*m|_1

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
RES0 = 0;
count =0;

% iterations
while(k <= params.Itnlim)
   
    % backtracking line-search
    t = t0;
    f0 = objFunc(m,dm,0,params);
    [f1, RES, TVterm] = objFunc(m,dm,t,params);
    lsiter = 0;
    while (f1 > f0 + alpha*t*(g0(:)'*dm(:))) && (lsiter<maxlsiter)
        t = t* beta;
        [f1, RES, TVterm] = objFunc(m,dm,t,params);
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
    fprintf('%d , relative residual: %f\n',...
            k, abs(RES-RES0)/RES);
    
    if (abs(RES-RES0)/RES <= gradToll); 
        count = count + 1;
    else
        count = 0;
    end 
    
    if (count == 10)
        break;
    end

    f = f1;
    RES0 = RES;
end

return;



function [obj,RES,TVterm] = objFunc(m, dm, t, params)
p = params.pNorm;
w1 = m+t*dm;

w2 = params.TV*w1;
TV = (w2.*conj(w2)+params.l1Smooth).^(p/2);
TVterm = sum(params.TVWeight(:).*TV(:));

RES = params.FT*w1 - params.data;
RES = RES(:)'*(params.wt(:).*RES(:));

obj = RES + TVterm;


function grad = wGradient(m,params)
p = params.pNorm;
w1 = params.TV*m;

gradTV = params.TV'*(p*w1.*(w1.*conj(w1)+params.l1Smooth).^(p/2-1));

gradRES = params.FT'*(params.wt.*((params.FT*m)-params.data));

grad = 2*gradRES + gradTV.*params.TVWeight;

