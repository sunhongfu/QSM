function[ extendedLocalPhase, sharpLocalPhase ] = esharp(totalPhase, mask, Options)
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
%   [eLP       = ESHARP(A,ROI) returns the 'extended local phase'
%   [eLP, sLP] = ESHARP(A,ROI) also returns 'sharp local phase'
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

numEcho  = size( totalPhase, 4 ) ;
numMasks = size( mask, 4 ) ;

if numMasks < numEcho
    disp('numMasks < numEchoes : using single mask for all echoes')
    mask = repmat( mask(:,:,:,1), [1 1 1 numEcho] ) ;
    
    disp('single set of Internal [expansion] Points (IP) will be used for extrapolation - FAST')
    isReusingIP = true ;

else
    isReusingIP = false ;
end


extendedLocalPhase = zeros( size( totalPhase ) ) ;
sharpLocalPhase    = zeros( size( totalPhase ) ) ;

if any( mod( size( totalPhase(:,:,:,1) ), 2 ) == 0)
    isCroppingToOddDimensions = true ;
    gridDimensionVector = size( totalPhase(:,:,:,1) ) - [1 1 1] ;
else
    isCroppingToOddDimensions = false ;
    gridDimensionVector = size( totalPhase(:,:,:,1) ) ;
end

%%sharp variables
midPoint            = ( gridDimensionVector + [1 1 1] ) / 2 ;
midPoint            = sub2ind( gridDimensionVector, midPoint(1), midPoint(2), midPoint(3) ) ;

sphere                  = createellipsoid( gridDimensionVector, Options.radii) ;
numAveragingPoints      = sum( sphere(:) ) ;

sharpFilter             = - sphere / numAveragingPoints ;
sharpFilter( midPoint ) = sharpFilter( midPoint ) + 1 ;

FFTSharpFilter          = fftc( sharpFilter ) ;


for echo = 1 : numEcho
    
    maskTE       = mask(:,:,:, echo) ;
    totalPhaseTE = maskTE .* totalPhase(:,:,:, echo) ;
    
    if isCroppingToOddDimensions
    totalPhaseTE   = makeodd( totalPhaseTE ) ;
    maskTE         = makeodd( maskTE ) ;
    end
    
    reducedROI   = shaver( maskTE, Options.radii ) ;
    
    %% SMV filtering
    tmpThresholdParameter   = 0.000001 ;
    
    % Deconv. (for eSHARP)
    tmp                   = fftc( reducedROI .* ifftc( fftc( sharpFilter ) .* fftc( totalPhaseTE) ) ) ./ FFTSharpFilter ;
    tmp( FFTSharpFilter < tmpThresholdParameter ) = 0 ;
    
    localPhaseBeforeSVD       = reducedROI .* real( ifftc( tmp ) );
    backgroundPhaseBeforeSVD  = reducedROI .* (totalPhaseTE - localPhaseBeforeSVD) ;
    
    
    % traditional SHARP
    tmp                   = fftc( reducedROI .* ifftc( fftc( sharpFilter ) .* fftc( totalPhaseTE) ) ) ./ FFTSharpFilter ;
    tmp( FFTSharpFilter < Options.thresholdParameter ) = 0 ;
    sharpLocalPhaseTE       = reducedROI .* real( ifftc( tmp ) );
    
    
    %% SHARP EDGES (Extrapolation)
    if (numEcho > 1) && (isReusingIP == true)
            if echo == 1
                Options.isSavingInterimVar = true ;
            else
                Options.isSavingInterimVar = false ;
                Options.IPtoEPDistance     = strcat(EdgeOut.Options.name,'_IPtoEPAssociations') ;
                Options.IPtoEdgeDistance   = strcat(EdgeOut.Options.name,'_IPtoEdgeDistance') ;
            end
    end
        
        
    EdgeOut  = sharpedges( backgroundPhaseBeforeSVD, maskTE, reducedROI, Options) ;
    
    extendedBackgroundPhase = EdgeOut.reducedBackgroundField + EdgeOut.extendedBackgroundField(:,:,:, EdgeOut.Options.expansionOrder + 1 ) ;
    
    tmpLocal       = totalPhaseTE - extendedBackgroundPhase ;
    fTmpLocal      = fftc( tmpLocal ) ;
    
    % Regularization
    fTmpLocal( FFTSharpFilter < Options.thresholdParameter )  = 0 ;
    
    localPhaseTSVD = maskTE .* ifftc( fTmpLocal) ;
    
    
    if isCroppingToOddDimensions
        extendedLocalPhase(:,:,:, echo) = makeodd(localPhaseTSVD, 'isUndoing') ;
        sharpLocalPhase(:,:,:, echo)    = makeodd(sharpLocalPhaseTE, 'isUndoing') ;
    else
        extendedLocalPhase(:,:,:, echo) = localPhaseTSVD ;
        sharpLocalPhase(:,:,:, echo)    = sharpLocalPhaseTE ;
    end
    
    
    
end
