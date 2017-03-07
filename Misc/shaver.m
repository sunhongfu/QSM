function[ ROIout ] = shaver( ROIin, R )
%SHAVER shaves a 3D binary mask by 'R'
%
%   SHAVER convolves binary input image ROI with an ellipsoid defined by
%   radius (or radii) R to return a contracted version of ROI
%
%   Syntax
%
%   ROIshaved = SHAVER(ROI,R)
%
%   Description
%
%   ROIshaved = SHAVER(ROI,R) returns ROI eroded by R. If R is a single
%   scalar, every dimension is eroded by R. If R is a 3-component vector
%   [Rx Ry Rz], each dimension is then eroded by its corresponding R value.
%

if any( R == 0 )
	ROIout = ROIin ;
	disp('Radius cannot be zero. Returning input array') 
else
	inputIsFloat = true;
    tmp = whos('ROIin');
	switch tmp.class
	    case {'single','double'}
	        inputIsFloat   = true ;
	        ROIin          = logical( ROIin ) ;
	    case {'logical'}
	        inputIsFloat = false ;
	end

	if any( mod( size( ROIin ), 2 ) == 0 )
	    isCroppingToOddDimensions = true ;
	    ROIin                     = makeodd( ROIin ) ;
	    gridDimensionVector       = size( ROIin ) ;
	else
	    isCroppingToOddDimensions = false ;
	    gridDimensionVector = size( ROIin ) ;
	end

	sphere              = createellipsoid( gridDimensionVector, R) ;
	ROIout              = ifftc( fftc(ROIin) .* fftc( sphere/sum( sphere(:)) ) ) ;
	ROIout              = abs(ROIout) >= 1 - 0.99/sum(sphere(:)) ;

	if inputIsFloat
	    ROIout = double( ROIout ) ;
	end

	if isCroppingToOddDimensions
	    ROIout = makeodd( ROIout, 'isUndoing' ) ;
	end

end
