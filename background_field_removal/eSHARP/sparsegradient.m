function[GradientTerms] = sparsegradient( dataArray, ROI, gradientOrder, voxelSize )
%SPARSEGRADIENT
%
%   SPARSEGRADIENT computes a succession of 3D gradients, storing desired
%   terms only, while minimizing the number of terms kept in memory at any
%   given moment without significantly compromising computation time
%
%   Syntax
%
%   SPARSEGRADIENT(A,ROI,order)
%   SPARSEGRADIENT(A,ROI,order,voxelSize)
%
%   Description
%
%   B = SPARSEGRADIENT(A,ROI,order) takes n = order derivatives of 3D array
%   A (default voxelSize = [1 1 1]), saving only terms at locations within
%   logial array ROI, and returns struct B with fields of
%
%       .coefficients 
%           gradient terms within ROI (each order (e.g. Gxxz) is a column
%           in a supersparse matrix)
%
%       .directions 
%           series of matrices definining gradient order and direction
%           (e.g. Gxxz is the 16th slice and looks like (1 0 0; 1 0 0; 0 0
%           1) )
%
%       .indexKey
%       only nonzero terms of B.coefficients (those belonging to ROI) are
%       stored: because the spatial coordinates of each of these cells is
%       the same for all orders, a 'sparse' matrix would be inefficient.
%       Instead, a single (sparse column vector) indexKey is returned: only
%       points belonging to ROI are nonzero - each containing the linear
%       index of the corresponding point in the compressed B.coefficients.
%       (e.g. an ROI consisting of two points P1>P2 would have
%       indexKey(P1)=1, indexKey(P2)=2)    
%
%           
% 2013
% topfer@ualberta.ca


DEFAULT_VOXELSIZE = [1 1 1] ;

%% check inputs

if nargin < 3 || isempty(dataArray) || isempty(ROI) || isempty(gradientOrder)
    error('Function requires at least three input arguments.')
end

if nargin < 4 || isempty(voxelSize)
    voxelSize = DEFAULT_VOXELSIZE ;
end
    


%% preliminaries
gridDimensionVector    = size( dataArray ) ;
ROI                    = logical( ROI(:) ) ;
numVoxels              = sum( ROI ) ;

% make indexKey...
GradientTerms.indexKey                = zeros( 1, prod(gridDimensionVector) ) ;
GradientTerms.indexKey( find( ROI ) ) = 1 : numVoxels ;
GradientTerms.indexKey                = sparse( GradientTerms.indexKey ) ;

for n = 0 : gradientOrder
    tmp(n + 1) = 3^n ;
end

gradInd(1, 1) = 1 ;

for n = 1 : gradientOrder + 1
    gradInd(n, 2)     = sum( tmp( 1:n ) ) ;
    gradInd(n + 1, 1) = gradInd(n, 2) + 1 ;
end
gradInd(gradientOrder + 2, :) = [] ;

numGradientTermsTotal = gradInd(gradientOrder + 1, 2);

% sparsified linear gradient terms (each term is a row in a sparse matrix):
GradientTerms.coefficients       = zeros( numVoxels, numGradientTermsTotal ); 
GradientTerms.coefficients(:, 1) = dataArray(ROI) ; 

GradientTerms.directions   = zeros( gradientOrder, 3, numGradientTermsTotal );

% the only 'full' gradient images stored @ a given time:
tmpCoeffs          = zeros([ gridDimensionVector gradientOrder + 1 ]) ; 
tmpCoeffs(:,:,:,1) = dataArray ;
currentTopography  = zeros( gridDimensionVector ) ;

numFrwdSteps  = gradientOrder - 1 ;
numBkwdSteps  = 0 ;
altitude      = 0 ;
lastStep      = 1 ;
pathTaken     = 0 ;
stepDirection = [1 0 0] ;

tmpDirections                = zeros(gradientOrder, 3) ;
n = 0 ;

if gradientOrder > 0


while lastStep > 0
    
    % Venture to top plateau: climbing along 'x'...
    for lastStep = 1 : numFrwdSteps
        
        % chart path
        altitude      = altitude + 1 ;
        lastStep      = altitude ;
        stepDirection = [1 0 0] ; % (i.e. 'x')
        tmpDirections( altitude, : ) = stepDirection ;
        stepDirectionStr       = char( stepDirection .* ['x' 'y' 'z'] ) ;
        route                  = find( stepDirectionStr ) ;
        pathTaken( altitude )  = route - 1 ; % (e.g. 0, 2, 1 = x -> z -> y)
        
        % measurement
        tmpCoeffs(:,:,:, lastStep + 1 ) = npointgradient( tmpCoeffs(:,:,:, lastStep), stepDirectionStr(route)) / voxelSize(route) ;
        
        % compression
        currentTopography          = tmpCoeffs(:,:,:, lastStep + 1) ;
        currentTopography          = currentTopography(ROI) ; 

        % GPS: translate coordinates
        lowerPlateau               = gradInd( altitude + 1, 1 ) ;
        currentLocation            = lowerPlateau + sum(pathTaken .* (3 .^( altitude  - 1: - 1 : 0) ));
        
        % Save measurement
        GradientTerms.directions(:,:, currentLocation)  = tmpDirections ;
        GradientTerms.coefficients(:, currentLocation ) = currentTopography ;
        
        % measurement index
        n = n + 1 ;
        
    end
    
    % check altimeter
    altitude     = altitude + 1 ;
    lastStep     = altitude ;
    lowerPlateau = gradInd( altitude + 1, 1 ) ;
    
    % Trident & Terminus: differentiate along -x,-y,-z, -> retreat
    for nextDirection = 1 : 3
        
        % chart path
        stepDirection                = [ 0 0 0 ] ;
        stepDirection(nextDirection) = 1 ;
        tmpDirections( altitude, : ) = stepDirection ;
        stepDirectionStr             = char( stepDirection .* ['x' 'y' 'z'] ) ;
        route                        = find( stepDirectionStr ) ;
        pathTaken(altitude)          = route - 1 ; % (e.g. 0, 2, 1 = x -> z -> y)
        
        % measurement
        tmpCoeffs(:,:,:, lastStep + 1)  = npointgradient( tmpCoeffs(:,:,:, lastStep), stepDirectionStr(route)) / voxelSize(route)  ;
        
        % compression
        currentTopography               = tmpCoeffs(:,:,:, lastStep + 1) ;
        currentTopography               = currentTopography(ROI) ; 
        
        % GPS: translate coordinates
        currentLocation            = lowerPlateau + sum(pathTaken .* (3 .^( altitude -1 : - 1 : 0) ));
        
        % Save measurement
        GradientTerms.directions(:,:, currentLocation)  = tmpDirections ;
        GradientTerms.coefficients(:, currentLocation ) = currentTopography ;
        
        % measurement index
        n = n + 1 ;
        
    end
    
    
    % Treacherous 'long -z: the road last traveled by (see: Frost, Robert. 1920)
    % to retreat: count consecutive passes along -z beginning w/latest
    numBkwdSteps = flipud(tmpDirections(:, 3))  ;
    numBkwdSteps = max(cumprod(numBkwdSteps) .* cumsum(numBkwdSteps)) ;
    
    % Retreat to camp...
    altitude                     = altitude - numBkwdSteps ;
    lastStep                     = altitude ;
    
    if altitude > 0
        
        % chart beginning of next path
        nextDirection                = find( tmpDirections(altitude, :) ) + 1 ;
        stepDirection                = [ 0 0 0 ] ;
        stepDirection(nextDirection) = 1 ;
        tmpDirections(altitude, :)   = stepDirection ;
        tmpDirections(altitude + 1 : gradientOrder, : ) = 0 ;
        
        % erase out-of-date part of chart
        pathTaken(altitude  : gradientOrder) = [] ;
        
        % update chart
        stepDirectionStr                = char( stepDirection .* ['x' 'y' 'z'] ) ;
        route                           = find( stepDirectionStr ) ;
        pathTaken(altitude)             = route - 1 ; % (e.g. 0, 2, 1 = x -> z -> y)
        
        % initial measurement of new set
        tmpCoeffs(:,:,:, lastStep + 1)  = npointgradient( tmpCoeffs(:,:,:, lastStep), stepDirectionStr(route)) / voxelSize(route) ;
        
        % compression
        currentTopography               = tmpCoeffs(:,:,:, lastStep + 1) ;
        currentTopography               = currentTopography(ROI) ; 
        
        % GPS
        lowerPlateau                    = gradInd( altitude + 1, 1 ) ;
        currentLocation                 = lowerPlateau + sum(pathTaken .* (3 .^( altitude - 1 : - 1 : 0) ));
        
        % save measurement
        GradientTerms.directions(:,:, currentLocation) = tmpDirections ;
        GradientTerms.coefficients(:,currentLocation ) = currentTopography ;
        
        % measurement index
        n = n + 1 ;
        
        % anticipate # steps to highest plateau (not peak)
        numFrwdSteps = numBkwdSteps -1 ;
        
    end
     
end

end


end

