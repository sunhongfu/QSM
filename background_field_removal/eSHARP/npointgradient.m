function pattern = npointgradient(pattern,direction,n)
%NPOINTGRADIENT fast calculation of gradient
%
%   NPOINTGRADIENTreturns gradient
%
%   Syntax
%
%   G = NPOINTGRADIENT(A,d,N)
%
%
%   Description
%   Returns gradient array G calculated by N-point central
%   differences of array A along direction d ('x' = 1st dim (i.e. row) ,
%   'y' = 2nd (i.e. column), 'z' = 3rd (i.e. slice) ).
%
%   Optional inputs
%       N (3, 5, 7, or 9)
%           default = 3
%
%   Easy derivation of formulas used @
%   http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/


DEFAULT_N = 3 ;

%% check input

if nargin < 2 || isempty(pattern) || isempty(direction)
    error('Function requires at least two input arguments.')
end

if nargin < 3 || isempty(n)
    n = DEFAULT_N ;
end

%%

if n == 3
    
    switch direction
        case 'x'
            pattern(2:end-1,:,:) = pattern(3:end,:,:) - pattern(1:end-2,:,:);
            pattern([1,end],:,:) = 0;
        case 'y'
            pattern(:,2:end-1,:) =pattern(:,3:end,:) - pattern(:,1:end-2,:) ;
            pattern(:,[1,end],:) = 0;
        case 'z'
            pattern(:,:,2:end-1) = pattern(:,:,3:end) - pattern(:,:,1:end-2) ;
            pattern(:,:,[1,end]) = 0;
    end
    pattern = pattern / 2 ;
    
elseif n == 5
    
    switch direction
        case 'x'
            pattern(3:end-2,:,:) = 8*(pattern(4:end-1,:,:) - pattern(2:end-3,:,:)) + pattern(1:end-4,:,:) - pattern(5:end,:,:) ;
            pattern([1,2,end-1,end],:,:) = 0;
        case 'y'
            pattern(:,3:end-2,:) = 8*(pattern(:,4:end-1,:) - pattern(:,2:end-3,:)) + pattern(:,1:end-4,:) - pattern(:,5:end,:) ;
            pattern(:,[1,2,end-1,end],:) = 0;
        case 'z'
            pattern(:,:,3:end-2) = 8*(pattern(:,:,4:end-1) - pattern(:,:,2:end-3)) + pattern(:,:,1:end-4) - pattern(:,:,5:end) ;
            pattern(:,:,[1,2,end-1,end]) = 0;
    end
    pattern = pattern / 12 ;
    
elseif n == 7
    switch direction
        case 'x'
            pattern(4:end-3,:,:) = 45*(pattern(5:end-2,:,:) - pattern(3:end-4,:,:)) + 9*(pattern(2:end-1,:,:) - pattern(6:end-1,:,:)) + pattern(7:end,:,:) - pattern(1:end-6,:,:)  ;
            pattern([1:3,end-2:end],:,:) = 0;
        case 'y'
            pattern(:,4:end-3,:) = 45*(pattern(:,5:end-2,:) - pattern(:,3:end-4,:)) + 9*(pattern(:,2:end-1,:) - pattern(:,6:end-1,:)) + pattern(:,7:end,:) - pattern(:,1:end-6,:)  ;
            pattern(:,[1:3,end-2:end],:) = 0;
        case 'z'
            pattern(:,:,4:end-3) = 45*(pattern(:,:,5:end-2) - pattern(:,:,3:end-4)) + 9*(pattern(:,:,2:end-1) - pattern(:,:,6:end-1)) + pattern(:,:,7:end) - pattern(:,:,1:end-6)  ;
            pattern(:,:,[1:3,end-2:end]) = 0;
    end
    pattern = pattern / 60 ;
    
elseif n == 9
    
    switch direction
        case 'x'
            pattern(5:end-4,:,:) = 672*(pattern(6:end-3,:,:) - pattern(4:end-5,:,:)) + 168*(pattern(7:end-2,:,:) - pattern(3:end-6,:,:)) + 32*(pattern(8:end-1,:,:) - pattern(2:end-7,:,:)) + 3*(pattern(9:end,:,:) - pattern(1:end-8,:,:)) ;
            pattern([1:4,end-3:end],:,:) = 0;
        case 'y'
            pattern(:,5:end-4,:) = 672*(pattern(:,6:end-3,:) - pattern(:,4:end-5,:)) + 168*(pattern(:,7:end-2,:) - pattern(:,3:end-6,:)) + 32*(pattern(:,8:end-1,:) - pattern(:,2:end-7,:)) + 3*(pattern(:,9:end,:) - pattern(:,1:end-8,:)) ;
            pattern(:,[1:4,end-3:end],:) = 0;
        case 'z'
            pattern(:,:,5:end-4) = 672*(pattern(:,:,6:end-3) - pattern(:,:,4:end-5)) + 168*(pattern(:,:,7:end-2) - pattern(:,:,3:end-6)) + 32*(pattern(:,:,8:end-1) - pattern(:,:,2:end-7)) + 3*(pattern(:,:,9:end) - pattern(:,:,1:end-8)) ;
            pattern(:,:,[1:4,end-3:end]) = 0;
    end
    pattern = pattern / 840 ;
    
end
end
