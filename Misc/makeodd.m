function[dataArray] = makeodd(dataArray, isPadding)

%MAKEODD returns dataArray with odd dimensions
%   
%   Syntax
%
%   A = MAKEODD(A)
%   A = MAKEODD(A, isPadding)
%
%   returns array A with odd dimensions
%   
%   Optional
%
%   isPadding
%       if false (DEFAULT), even dimensions have a slice shaved off
%       if true, even dimensions are padded with a slice of zeros
%

DEFAULT_ISPADDING = false ;

if nargin < 2
    isPadding = DEFAULT_ISPADDING ;
end


gridDimensionVector = size( dataArray ) ;
isOdd               = mod( gridDimensionVector, 2 ) ;
isEven              = double(~isOdd) ;


if isPadding
    
    dataArray = padarray(dataArray, isEven, 'post') ;
    
else
   
    oddLimits = gridDimensionVector - isEven ;
if numel( oddLimits) == 1
	dataArray = dataArray( 1:oddLimits(1) ) ;
elseif numel( oddLimits) ==2
	dataArray = dataArray( 1:oddLimits(1), 1:oddLimits(2) ) ;
else
    dataArray = dataArray( 1:oddLimits(1), 1:oddLimits(2), 1:oddLimits(3) ) ;
    
end


end
        




