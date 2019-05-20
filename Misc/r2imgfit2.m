function [R2 T2 amp] = r2imgfit2(img,te,Wt)
%   Function to fit exponentials to an image
%   
%   Author: Marc Lebel 08/2009
%   Usage:
%   [R2 T2 amp off] = r2fit_wrap(img,TE,Wt)
%   Input:
%   img = image of size M x N x slices x nTE
%   TE = array of echo times (1 x nTE) (s)
%   Wt =  a weighting function of size M x N x slices x nTE
%   
%   Output:
%   R2   = Decay rates (1/s)
%   T2   = Decay times (s)
%   amp  = amplitudes

%   Check inputs
if nargin < 2
    error('Function needs at least 2 inputs');
end
if nargin < 3 || isempty(Wt)
    Wt = ones(size(img));
end

%   Get image size
[np nv ns ne] = size(img);
if ne ~= length(te)
    error('Inconsistent data');
end
if ~isreal(img)
    img = abs(img);
end

%   Initialize output
R2  = zeros(np,nv,ns);
T2  = zeros(np,nv,ns);
amp = zeros(np,nv,ns);

%   Normalize image
[a,b] = hist(img(:),1024);
th = b(a == max(a));clear a b
img = img - th;
img = img./mean(img(:));

%   Define fit options
opt = optimset('MaxIter',50,'MaxFunEvals',200,'TolFun',1e-4,...
    'TolCon',0,'TolX',1e-3,'LargeScale','off','Display','off',...
    'DiffMaxChange',1e1,'DiffMinChange',1e-5,'Algorithm','active-set',...
    'FunValCheck','off');

%   Loop through image
poolobj = parpool;

parfor k=1:ns
%for k = 43:43
    ampt = zeros(np,nv);
    R2t = zeros(np,nv);
    T2t = zeros(np,nv);
    imgt = img(:,:,k,:);
    Wtt = Wt(:,:,k,:);
    for i=1:np
    for j=1:nv
        S = imgt(i,j,1,:);
        W = Wtt(i,j,1,:);
        if any(S>0) && any(W>0) && all(S(1:2)>0.5)
            
            %   Define fit bounds
            %    [   A    T2 ];
            X0 = [ S(1) 0.075];
            Xl = [ 0.0  0.005];
            Xu = [20.0  4.000];
            
            X = fmincon(@wexpval,X0,[],[],[],[],Xl,Xu,[],opt,S(:),te(:),W(:));
            ampt(i,j) = X(1);
            R2t(i,j) = 1./(X(2)+eps);
            T2t(i,j) = X(2);
            %plot(te,S(:),'o',te,X(1)*exp(-te./X(2)));drawnow;
        end
    end
    end
    amp(:,:,k) = ampt;
    R2(:,:,k) = R2t;
    T2(:,:,k) = T2t;
    fprintf('Done slice %g of %g\n',k,ns);
end

delete(poolobj)

return

function SSE = wexpval(X,S,TE,W)
SSE = sum( (W .* ((X(1)*exp(-TE./X(2))) - S)).^2);
%plot(TE,S,'o',TE,(X(1)*exp(-TE./X(2))));drawnow;
return
