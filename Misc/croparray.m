function[ dataArray ] = croparray( dataArray, outputSize ) 

%CROPARRAY returns cropped/central portion of 3D array
%
% Syntax
%   B = croparray( A, desiredSize )
%   
%   returns central portion of 3D array A: data outside 'desiredSize' is
%   trimmed away.
%
% See Also
% padarray

% TODO
% should work for arbitrary dimensions
% odd & even arrays
% check outputSize < inputSize

gridDimensionVector = size( dataArray ) ;


isOdd               = logical(mod( gridDimensionVector, 2 )) ;

if any(isOdd)
% if odd
midPoint = (gridDimensionVector + [1 1 1])/ 2 ; 

low = midPoint - outputSize/2 + [1 1 1] ;
high = midPoint + outputSize/2 ;

dataArray = dataArray( low(1):high(1), low(2):high(2), low(3):high(3)) ;

else

low = gridDimensionVector/2 - outputSize/2 + [1 1 1];
high = gridDimensionVector/2 + outputSize/2 ;

dataArray = dataArray( low(1):high(1), low(2):high(2), low(3):high(3) ) ;



end

