function[ unitDipoleField ] = createunitdipole( gridDimensionVector, Options )
%CREATEUNITDIPOLE creates unit dipole kernel
%
%   CREATEUNITDIPOLE  
%
%   Syntax
%
%   CREATEUNITDIPOLE( gridDimensionVector ) 
%   CREATEUNITDIPOLE( gridDimensionVector, Options )
%
%   Description
%
%   unitDipoleField = CREATEUNITDIPOLE( gridDimensionVector, Options )
%   returns unit dipole field/kernel of size given by gridDimensionVector
%
%   The following Option-fields are supported
%       .voxelSize 
%           k-space step size becomes 
%           dK = (Options.gridDimensionVector .* voxelSize)^-1 
%           default = [1 1 1]

% TODO
% spatial or Fourier generating options
%
% FIX odd/even dimensions. should be able to cope w/odd-x even-y for
% example.
% This function should call 'makekspace' to form coordinate space


%%
DEFAULT_RETURNDOMAIN = 'Spatial' ;
DEFAULT_VOXELSIZE   = [ 1 1 1] ;

Options.returnDomain = DEFAULT_RETURNDOMAIN ;

if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end

FOV  = gridDimensionVector .* Options.voxelSize ;
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

%%

isOdd               = logical(mod( gridDimensionVector, 2 )) ;

if isOdd
    
    midPoint        = (gridDimensionVector + [1 1 1])/2 ;
    
    limits          = midPoint - [1 1 1] ;
    
    [kx, ky, kz]    = ndgrid( (-limits(1):limits(1))/FOVx, (-limits(2):limits(2))/FOVy, (-limits(3):limits(3))/FOVy ) ;
    
    unitDipoleField = 1/3 - (kz .^2 )./ ( kx .^2 + ky .^2 + kz .^2 );
        
    unitDipoleField( midPoint(1), midPoint(2), midPoint(3) ) = 0 ;
    
    unitDipoleField = ifftshift( unitDipoleField ) ;
    
    switch Options.returnDomain
        case 'Spatial'
            unitDipoleField = ifftc( unitDipoleField ) ;
    end
    
else
    
    limits = gridDimensionVector/2 ;
        
    [kx,ky,kz] = ndgrid((-limits(1):limits(1)-1)/FOVx, (-limits(2):limits(2)-1)/FOVy, (-limits(3):limits(3)-1)/FOVz) ;
    
    unitDipoleField = 1/3 - (kz .^2 )./(kx.^2 + ky.^2 + kz.^2);
    
    unitDipoleField = ifftshift(unitDipoleField);
    
    unitDipoleField(1,1,1) = 0;
    
    switch Options.returnDomain
        case 'Spatial'
            unitDipoleField = ifftc( unitDipoleField) ;
            
    end
end

end




