function[ ellipsoid ] = createellipsoid(gridDimensions, radii, offset)
%CREATEELLIPSOID 
%
%   Creates ellipsoidal region (cells == 1) within an otherwise 0 array
%
%
%   Syntax
%
%   CREATEELLIPSOID(gridDimensions, radii)
%   CREATEELLIPSOID(gridDimensions, radii, offset)
%
%
%   Description
%
%   ellipsoid = CREATEELLIPSOID(gridDimensions, radii) returns an
%   array of size gridDimensions with an ellipsoid defined by radii at
%   the center
%
%   ellipsoid =CREATEELLIPSOID(gridDimensions, radii, offset) center
%   of the ellipsoid is offset
% 
%
%   all inputs may be single scalars or 3 element vectors (a scalar
%   will treat every direction the same: scalarValue *[ 1 1 1 ] )
%   
%   gridDimensions should consist of odd numbers
%
%

% check arguments
if nargin < 3
    offset = [0 0 0];
end

if length(gridDimensions) == 1
    gridDimensions = gridDimensions * [ 1 1 1 ];
end

if length(radii) == 1
    radii = radii * [ 1 1 1 ] ;
end

if length(offset) == 1
    offset = offset * [ 1 1 1 ] ;
end



% creates ellipsoid
[x,y,z]   = ndgrid(-radii(1) : radii(1), -radii(2) : radii(2), -radii(3) : radii(3)) ; % coordinates
ellipsoid = ( x .^2 /radii(1)^2 + y .^2 /radii(2)^2 + z .^2 /radii(3)^2 ) <= 1 ; % equation of ellipsoid

% places it in larger array
temp  = zeros( gridDimensions );
origin = ( gridDimensions + [1 1 1] ) /2 ;
offset = origin + offset;
temp( offset(1) - radii(1) : offset(1) + radii(1), offset(2) - radii(2) : offset(2) + radii(2), offset(3) - radii(3) : offset(3) + radii(3)) = ellipsoid;

ellipsoid = temp;

end