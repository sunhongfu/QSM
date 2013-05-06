function[ localPhase ] = esharp(totalPhase, mask, Options)
%ESHARP - extended SHARP
%
%   ESHARP returns the local field map over the ROI provided
%
%
%   Syntax
%
%   ESHARP(A,ROI)
%   ESHARP(A,ROI,Options)
%
%
%   Description
%
%   B = ESHARP(A,ROI) returns 
%   
%   [B, outputName] = ESHARP(A,ROI,reducedROI)
%       also returns the filename where B was written to disk (done by
%       default - to not bother writing to disk, set
%       Options.isWritingToDisk = 'false' )
%   
%   The following Option-fields are supported
%
%       voxelSize
%               default: [1 1 1] (isotropic)
%
%       thresholdParameter
%           SHARP regularization (TSVD) parameter
%               default: 0.05
%
%
%       radii
%           radii of SMV kernel
%               default: [6 6 6]./voxelSize
%               
%
%   See HELP SHARPEDGES for other Option-fields
%
%   R Topfer 2012     topfer@ualberta.ca



%% constants
DEFAULT_VOXELSIZE                       = [1 1 1] ;
DEFAULT_RADII                           = [6 6 6] ; 
DEFAULT_THRESHOLDPARAMETER              = 0.05 ;



%% check inputs

if nargin < 2 || isempty(totalPhase) || isempty(mask) 
    error('Function requires at least 2 input arguments.')
end

if nargin < 3 || isempty(Options)
    disp('Default parameters will be used')
    disp('assuming isotropic resolution')
    Options.dummy = [] ;
end



if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end

if ~myisfield( Options, 'radii' ) || isempty(Options.radii)
    Options.radii = round( DEFAULT_RADII ./ Options.voxelSize )  ;
end
  
if ~myisfield( Options, 'thresholdParameter' ) || isempty(Options.thresholdParameter)
    Options.thresholdParameter = DEFAULT_THRESHOLDPARAMETER ;
end


%% Initialize


if any( mod( size( totalPhase ), 2 ) == 0)
    disp('requires data array with odd-dimensions: cropping data')
    totalPhase = makeodd( totalPhase ) ;
    mask       = makeodd( mask ) ;
end

totalPhase = - mask .* totalPhase;

gridDimensionVector = size( totalPhase ) ;
midPoint            = ( gridDimensionVector + [1 1 1] ) / 2 ;
midPoint            = sub2ind( gridDimensionVector, midPoint(1), midPoint(2), midPoint(3) ) ;



%% SMV filtering
tmpThresholdParameter   = 0.000001 ; 
reducedROI              = shaver( mask, Options.radii ) ;

sphere                  = createellipsoid( gridDimensionVector, Options.radii) ;
numAveragingPoints      = sum( sphere(:) ) ;

sharpFilter             = - sphere / numAveragingPoints ;
sharpFilter( midPoint ) = sharpFilter( midPoint ) + 1 ;

FFTSharpFilter          = fftc( sharpFilter ) ;

% Deconv.
tmp                        = fftc( reducedROI .* ifftc( fftc( sharpFilter ) .* fftc( totalPhase) ) ) ./ FFTSharpFilter ;
tmp( FFTSharpFilter < tmpThresholdParameter ) = 0 ;

localPhaseNoSVD       = reducedROI .* real( ifftc( tmp ) );
backgroundPhaseNoSVD  = reducedROI .* (totalPhase - localPhaseNoSVD) ;



%% SHARP EDGES (Extrapolation)
EdgeOut  = sharpedges( backgroundPhaseNoSVD, mask, reducedROI, Options) ;

extendedBackgroundPhase = EdgeOut.reducedBackgroundField + EdgeOut.extendedBackgroundField(:,:,:, EdgeOut.Options.expansionOrder + 1 ) ;

tmpLocal       = totalPhase - extendedBackgroundPhase ;
fTmpLocal      = fftc( tmpLocal ) ;

% Regularization
fTmpLocal( FFTSharpFilter < Options.thresholdParameter )  = 0 ;

localPhaseTSVD = mask .* ifftc( fTmpLocal) ;


localPhase = localPhaseTSVD ;




