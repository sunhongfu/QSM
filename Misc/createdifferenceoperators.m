function [Dx, Dy, Dz] = createdifferenceoperators( gridSize, gridSpacing, nOrder )
%CREATEDIFFERENCEOPERATORS
%
%   Returns sparse ("central difference") matrices Dx, Dy, & Dz that act as partial
%   differential operators in row, column, and slice directions respectively.
%
%   Syntax
%
%   [Dx, Dy, Dz] = CREATEDIFFERENCEOPERATORS( gridSize, gridSpacing, order )
%   
%   Input
%   
%   gridSize 
%       Actual 3D size of the vector the differential operators will 
%       ultimately be applied to (i.e. gridSize = [nRows nColumns nSlices] )
%
%   gridSpacing
%       distance b/tw lattice points [dx dy dz]
%
%   order 
%       == 1
%           1st order fwd differences 
%           e.g. Dx*b + Dy*b + Dz*b = gradient of b
%
%       == 2 
%           2nd order central differences 
%           e.g. Dx*b + Dy*b + Dz*b = Laplacian of b
%   
%    
% 2014
% topfer@ualberta.ca

%% check inputs
errorMsg = 'Check usage: type HELP CREATEDIFFERENCEOPERATORS' ;
assert( nargin == 3, errorMsg ) ; 
assert( ((nOrder > 0) && (nOrder < 3)), errorMsg ) ;
assert( numel(gridSize) == 3, errorMsg );
assert( all(gridSize > 0), errorMsg );
assert( all(isfinite(gridSize) ), errorMsg ) ;


%% initialize
N = prod( gridSize ) ;

n1 = gridSize(1) ;
n2 = gridSize(2) ;
n3 = gridSize(3) ;

numCellsPerSlice = n2 * n1 ;

if nOrder == 1
%%
    % 1st dim
    e  = ones( n1, 1 ) ;
    
    Dx = spdiags([-e e], [-1 1], n1, n1) ;
    while size(Dx,2) < N
        Dx = blkdiag(Dx, Dx) ;
    end
    
    Dx = Dx(1 : N, 1 : N) ;
    
    % 2nd dim
    e = ones( numCellsPerSlice, 1 ) ;
    
    Dy = spdiags([-e e], [-n1 n1], numCellsPerSlice, numCellsPerSlice) ;
    while size(Dy,2) < N
        Dy = blkdiag(Dy, Dy) ;
    end
    
    Dy = Dy(1 : N, 1 : N) ;
    
    
    % 3rd dim
    e = ones( N, 1) ;
    Dz = spdiags([-e e], [-numCellsPerSlice numCellsPerSlice], N, N) ;
    
    Dx = Dx/(2*gridSpacing(1)) ;
    Dy = Dy/(2*gridSpacing(2)) ;
    Dz = Dz/(2*gridSpacing(3)) ;
   

elseif nOrder == 2
    
    % 1st dim
    e  = ones( n1, 1 ) ;

    Dx = spdiags([e -2*e e], -1:1, n1, n1) ;
    while size(Dx,2) < N
        Dx = blkdiag(Dx, Dx) ;
    end
    
    Dx = Dx(1 : N, 1 : N) ;
    
    % 2nd dim
    e = ones( numCellsPerSlice, 1 ) ;

    Dy = spdiags([ e -2*e e], [-n1 0 n1], numCellsPerSlice, numCellsPerSlice) ;
    while size(Dy,2) < N
        Dy = blkdiag(Dy, Dy) ;
    end
    
    Dy = Dy(1 : N, 1 : N) ;
    
    % 3rd dim
    e = ones( N, 1) ;
    Dz = spdiags([e -2*e e], [-numCellsPerSlice 0 numCellsPerSlice], N, N) ;

    Dx = Dx/(gridSpacing(1)^2) ;
    Dy = Dy/(gridSpacing(2)^2) ;
    Dz = Dz/(gridSpacing(3)^2) ;

end
