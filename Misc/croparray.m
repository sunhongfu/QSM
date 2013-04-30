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

% if odd
midPoint = (gridDimensionVector + [1 1 1])/ 2 ; 

tmp1 = midPoint - (outputSize + [1 1 1])/2 + [1 1 1] ;
tmp2 = midPoint + (outputSize + [1 1 1])/2 - [1 1 1] ;

dataArray = dataArray( tmp1(1):tmp2(1), tmp1(2):tmp2(2), tmp1(3):tmp2(3)) ;

end

