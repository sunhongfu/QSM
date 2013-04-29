function[kArray]=makekspace(gridDimensionVector, Options) 
%MAKEKSPACE produces 3D array of k-space coordinates
%
%   MAKEKSPACE returns K(kx,ky,kz)  
%
%   Syntax
%
%   MAKEKSPACE(gridDimensionVector)
%   MAKEKSPACE(gridDimensionVector,Options)
%
%   Description
%
%   K = MAKEKSPACE(gridDimensionVector,Options) returns k-space coordinate
%   array
%
%   The following Option-fields are supported
%       .voxelSize 
%           k-space step size becomes 
%           dK = 2pi*(Options.gridDimensionVector .* voxelSize)^-1 
%           default = [1 1 1]

% TODO
% FIX odd/even dimensions. should be able to cope w/odd-x even-y for
% example.

%% constants

DEFAULT_VOXELSIZE                       = [1 1 1] ;


%% check inputs
if nargin < 2 || isempty(Options)
    disp('Assuming isotropic resolution...')
    Options.dummy = [] ;
end

if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end

%% intialize
kArray = zeros([gridDimensionVector 3]);
FOV    = gridDimensionVector .* Options.voxelSize ;
FOVx   = FOV(1)/(2*pi);
FOVy   = FOV(2)/(2*pi);
FOVz   = FOV(3)/(2*pi);


isOdd  = logical(mod( gridDimensionVector, 2 )) ;

if isOdd
    
    midPoint        = (gridDimensionVector + [1 1 1])/2 ;
    limits          = midPoint - [1 1 1] ;  
    [kx, ky, kz]    = ndgrid( (-limits(1):limits(1))/FOVx, (-limits(2):limits(2))/FOVy, (-limits(3):limits(3))/FOVy ) ;
    
else
    
    limits       = gridDimensionVector/2 ;
    [kx, ky, kz] = ndgrid((-limits(1):limits(1)-1)/FOVx, (-limits(2):limits(2)-1)/FOVy, (-limits(3):limits(3)-1)/FOVz) ;
    
end

kArray(:,:,:,1) = kx ;
kArray(:,:,:,2) = ky ;
kArray(:,:,:,3) = kz ;

end
