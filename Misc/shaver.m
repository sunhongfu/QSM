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
%   Ryan Topfer topfer@ualberta.ca
% gitTEST
sphere              = createellipsoid( size( ROIin ), R) ;
ROIout              = ifftc( fftc(ROIin) .* fftc( sphere/sum( sphere(:)) ) ) ;
ROIout              = abs(ROIout) >= 1 - 0.99/sum(sphere(:)) ;
end
