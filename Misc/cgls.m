function [x] = cgls(leastSquaresOperator, b, x, Options)
%CGLS Conjugate gradient for least-squares problems
%
% Ax = b
%
%......
%   TODO
%   
%       documentation
%
%       fix self-adjoint conditional

%% constants
DEFAULT_MAXITERATIONS = 500 ;
DEFAULT_TOLERANCE     = 1E-6 ;
DEFAULT_ISDISPLAYINGPROGRESS = false ;
%% check inputs

if nargin < 2 || isempty(leastSquaresOperator) || isempty(b)
    error('CGLS requires at least 2 arguments')
end

if nargin == 2 || isempty(x)
    x = 0;
end

if nargin < 4 || isempty( Options )
    Options.dummy = [] ;
end

if  ~myisfield( Options, 'maxIterations' ) || isempty(Options.maxIterations)
    Options.maxIterations = DEFAULT_MAXITERATIONS ;
end

if  ~myisfield( Options, 'tolerance' ) || isempty(Options.tolerance)
    Options.tolerance = DEFAULT_TOLERANCE ;
end

if  ~myisfield( Options, 'isDisplayingProgress' ) || isempty(Options.isDisplayingProgress)
    Options.isDisplayingProgress = DEFAULT_ISDISPLAYINGPROGRESS ;
end





Tmp = whos('leastSquaresOperator') ;

switch Tmp.class
    case 'function_handle'
        applyLSOperator = leastSquaresOperator ;
    
    case 'double'
        applyLSOperator = @applymatrixoperator ;
end

  
isSelfAdjoint = false ;    


%% solve
s           = b - applyLSOperator(x) ;
r           = applyLSOperator(b - applyLSOperator(x), 1) ;
p           = r ;
discrepancy = 1 ;
k           = 0 ;
normB       = norm( b(:) ) ;

while k < Options.maxIterations && discrepancy > Options.tolerance
    
    q        = applyLSOperator(p);
    
    rtr      = r(:)' * r(:) ;
    qtq      = q(:)' * q(:) ;
    alpha    = rtr / qtq ;
    x        = x + alpha*p ;
    s        = s - alpha*q ;
    
    r        = applyLSOperator(s,1) ;
    
    beta     = r(:)'*r(:)/(rtr) ;
    p        = r + beta*p ;
    
    Ax          = applyLSOperator(x) ;
    discrepancy = norm(Ax(:) - b(:))/normB ;
    
    if Options.isDisplayingProgress
        disp( [int2str(k) ' ' num2str(discrepancy) ] )
    end

    k = k + 1 ;
end
disp( ['Discrepancy after ' int2str(k) ' iterations: '  num2str(discrepancy) ] )
    clear matrixOperator


function b = applymatrixoperator(x, isAdjoint)
    
    if (nargin > 1) && (isAdjoint == 1) && (~isSelfAdjoint)  
        b = leastSquaresOperator' * x;
    else
        b = leastSquaresOperator * x ;
    end
end


        
end