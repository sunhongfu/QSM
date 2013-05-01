function[ EdgeOut ] = sharpedges(dataArray, ROI, reducedROI, Options)

%SHARPEDGES recovery of phase information near ROI boundary post-SHARP
%
%   SHARPEDGES returns the background field map over the entire ROI
%
%
%   Syntax
%
%   SHARPEDGES(A,ROI,reducedROI)
%   SHARPEDGES(A,ROI,reducedROI,Options)
%
%
%   Description
%
%   B = SHARPEDGES(A,ROI,reducedROI) returns object containing the
%   following variables:
%   
%    .......................
%   
%   The following Option-fields are supported
%       expansionOrder
%               default: 2
%
%       voxelSize
%               default: [1 1 1] (isotropic)
%
%       aspectRatio
%               default: determined by voxelSize 
%
%       name
%           (of .mat where variables are saved)
%               default: 'sharpEdgesResults'
%
%       numIterations
%           convergence is examined after the first iteration: points for
%           which the series didn't intially converge will be re-evaluated
%           at iteration 2 after incorporating the points for which the
%           solver did converge
%               default: 1
%
%       numCPU
%           specifies number of processing units used for parallel
%           computation.
%               default: matlabpool ~ whatever is available
%
%       isUsingConvergenceCriterion
%           points for which the series does not appear to be converging
%           after the 1st iteration will be solved for again in the 2nd
%           iteration using different IEPs. 
%               default: false
%       
%       The following options are more to save time during debugging than anything
%       else since they essentially require running the code once in the
%       1st place:
%
%       For all of the below, the default is simply to determine
%       them during run-time. 
%       
%
%       ALSO: Although these are essentially 3D things, they are SAVED in a
%       COMPRESSED LINEAR form
%
%       IEPtoEdgeDistance
%       specifies path of MATLAB matrix containing the variables
%           -IEPtoEdgeDistance : determines 'harmonic neighbourhood' about
%           each IEP
%               
%       IEPtoEPassociations
%           specifies path of MATLAB matrix containing the variables
%           -associatedIEPs : linear indices of IEPs for each EP 
%           
%           -associatedIEPtoEPDistance : Euclidean distance b/tw them
%           
%           -associatedIEPtoEdgeDistance : Each EP contains the
%           IEPtoEdgeDistance of its corresponding IEP.               
%
%       GradientTerms
%           specifies path of MATLAB matrix containing the variables
%           -GradientTerms : Object output from sparsegradient
%
%       taylorTerms
%           specifies path of MATLAB matrix containing the variables
%           -taylorTerms : expansionOutput
%
% ABBRV:
%       IEP = Internal Expansion Point : point at which extrapolation is based   
%       EP  = Extension Point : looking to determine Bkgr. Field here 
%
%
% R Topfer 2012     topfer@ualberta.ca
%


disp(' ')
ComputationTime.total = tic ;

%% constants

DEFAULT_NAME                            = 'sharpEdgesResults' ;
DEFAULT_EXPANSIONORDER                  = 2 ;
DEFAULT_NUMITERATIONS                   = 1 ;
DEFAULT_NUMCPU                          = findResource ;
DEFAULT_NUMCPU                          = DEFAULT_NUMCPU.ClusterSize ;
DEFAULT_VOXELSIZE                       = [1 1 1] ;

DEFAULT_ISUSINGCONVERGENCECRITERION     = false ;

DEFAULT_ISMAPPINGEXPANSIONPOINTS        = true ;
DEFAULT_ISMAPPINGHARMONICNEIGHBOURHOODS = true ;
DEFAULT_ISCALCULATINGGRADIENTS          = true ;
DEFAULT_ISPERFORMINGTAYLOREXPANSION     = true ;

DEFAULT_ISWRITINGTODISK                 = true ;

DEFAULT_ISDEBUGGING                     = false;
DEFAULT_ISSAVINGINTERIMVAR              = false ;

%% check inputs

if nargin < 3 || isempty(dataArray) || isempty(ROI) || isempty(reducedROI)
    error('Function requires at least three input arguments.')
end

if nargin < 4 || isempty(Options)
    disp('Default parameters will be used')
    Options.dummy = [] ;
end

if  ~myisfield( Options, 'name' ) || isempty(Options.name)
    Options.name = DEFAULT_NAME ;
end

if  ~myisfield( Options, 'expansionOrder' ) || isempty(Options.expansionOrder)
    Options.expansionOrder = DEFAULT_EXPANSIONORDER ;
end

if  ~myisfield( Options, 'offset' ) || isempty(Options.offset)
        Options.offset = Options.expansionOrder + 1 ; % How 'far back' from the edge of the eroded mask the expansion points are
%     Options.offset = Options.expansionOrder + 2 ; % How 'far back' from the edge of the eroded mask the expansion points are
end

if  ~myisfield( Options, 'numIterations' ) || isempty(Options.numIterations)
    Options.numIterations = DEFAULT_NUMITERATIONS ;
end

if  ~myisfield( Options, 'numCPU' ) || isempty(Options.numCPU)
    Options.numCPU = DEFAULT_NUMCPU ;
end

if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end

if ~myisfield( Options, 'aspectRatio' ) || isempty(Options.aspectRatio)
    Options.aspectRatio = round( Options.voxelSize / min( Options.voxelSize ) ) ;
end

if ~myisfield( Options, 'isUsingConvergenceCriterion') || isempty(Options.isUsingConvergenceCriterion)
    Options.isUsingConvergenceCriterion = DEFAULT_ISUSINGCONVERGENCECRITERION ;
end

if  ~myisfield( Options, 'IEPtoEPassociations' ) || isempty(Options.IEPtoEPassociations)
    isMappingExpansionPoints = DEFAULT_ISMAPPINGEXPANSIONPOINTS ;
else
    isMappingExpansionPoints = false ;
end

if ~myisfield( Options, 'IEPtoEdgeDistance') || isempty(Options.IEPtoEdgeDistance)
    isMappingHarmonicNeighbourhoods = DEFAULT_ISMAPPINGHARMONICNEIGHBOURHOODS ;
else
    isMappingHarmonicNeighbourhoods = false ;
end
    
if ~myisfield( Options, 'GradientTerms') || isempty(Options.GradientTerms)
    isCalculatingGradients = DEFAULT_ISCALCULATINGGRADIENTS ;
else
    isCalculatingGradients = false ;
end

if ~myisfield( Options, 'taylorTerms') || isempty(Options.taylorTerms)
    isPerformingTaylorExpansion = DEFAULT_ISPERFORMINGTAYLOREXPANSION ;
else
    isPerformingTaylorExpansion = false ;
end

if ~myisfield( Options, 'isDebugging') || isempty(Options.isDebugging)
    Options.isDebugging = DEFAULT_ISDEBUGGING ;
end

if ~myisfield( Options, 'isSavingInterimVar') || isempty(Options.isSavingInterimVar)
    Options.isSavingInterimVar = DEFAULT_ISSAVINGINTERIMVAR ;
end

if ~myisfield( Options, 'isWritingToDisk') || isempty(Options.isWritingToDisk)
    Options.isWritingToDisk = DEFAULT_ISWRITINGTODISK ;
end

    
Options.name = strcat( Options.name, '_order', int2str(Options.expansionOrder), '_iteration', int2str(0) );


%% Common to all iterations...
% constant terms
gridDimensionVector      = size( dataArray ) ;



%% CHANGE THIS 
if sum( mod(gridDimensionVector,2) ) < 3
    error('input dimensions must be ODD') ;
end




%%

inputField               = dataArray ;

ROIedge                  = ROI - shaver(ROI, 1) ;
edgeIndices              = find( ROIedge ) ;
numEdgePoints            = length( edgeIndices ) ;
[X1edge, Y1edge, Z1edge] = ind2sub( gridDimensionVector, edgeIndices ) ;




%% For the first iteration...
reducedBackgroundField  = reducedROI .* inputField ;
previousIEPs            = zeros( gridDimensionVector ) ;

% The following determines the number of gradient terms required for
% each order (essentially the 2nd column of gradientIndices)
% The 1st column demarcates the index at which that particular order of
% gradient begins (e.g Fo, Fx, Fy, Fz, Fxx, Fxy, etc. - order 1 begins at
% index 2 and ends at index 4, order 2 begins at index 5 and ends at index
% 13, etc.)

for n = 0 : Options.expansionOrder
    tmp(n + 1) = 3^n ;
end

gradInd(1, 1) = 1 ;

for n = 1 : Options.expansionOrder + 1
    gradInd(n, 2)     = sum( tmp( 1:n ) ) ;
    gradInd(n + 1, 1) = gradInd(n, 2) + 1 ;
end

gradInd(Options.expansionOrder + 2, :) = [] ;

numGradientTermsTotal = gradInd(Options.expansionOrder + 1, 2) ;
    


%% Iterative stuff
for currentIteration = 1 : Options.numIterations
    
    % Reset
    GradientTerms           = [] ;
    convergenceSurface      = [] ;
    expansionCapacity       = [] ;
    extendedBackgroundField = [] ;
    
        
    % Update output file name
    previousIteration = currentIteration - 1 ;
    pIStr             = int2str( previousIteration ) ;
    cIStr             = int2str( currentIteration ) ;
    lName             = length(Options.name) ;
    lPrevStr          = length(pIStr) ;
    Options.name( lName - lPrevStr + 1 : lName - lPrevStr + length(cIStr) ) = cIStr ;
    
    % update extendedROI
    extendedROI          = ROI - reducedROI ;
    extensionPoints      = find( extendedROI ) ;
    numPointsExtendedROI = length( extensionPoints ) ; % num of pt.s within extendedROI
    [X1, Y1, Z1]         = ind2sub( gridDimensionVector, extensionPoints ) ;
    
    % update expansionSurface
    tmpROI           = reducedROI - reducedROI .* ROIedge ;  % for iterations >1, unsure what the implications of including edges might be
    tmpROI           = shaver(tmpROI, Options.offset ) ; % back up from edges by offset...
    tmpROIedge       = tmpROI - shaver(tmpROI, 1 ) ; % avoid using the same IEPs one iteration to the next...                      
    tmpROI           = tmpROI - tmpROIedge .* previousIEPs ; % concern is exclusively w/previous IEPs on tmpROIedge : w/out this, in the next line, we'd introduce more pt.s into the pool of potential IEPs that would not be chosen/used but their inclusion would slow down mapping steps  
    tmpROI           = tmpROI - shaver(tmpROI, 1 ) ; % reduced to ~ lamina
    expansionSurface = tmpROI(:) ;
    
    clear tmpROI tmpROIedge
    
    internalExpansionPoints = find( expansionSurface ) ;
    numIEPs                 = length( internalExpansionPoints ) ; % num of pt.s to choose from for determining best expansion point
    [Xo, Yo, Zo]            = ind2sub( gridDimensionVector, internalExpansionPoints ) ;
    
        
    disp( [ 'Solving for extended background field: Iteration ' int2str(currentIteration) ])
    disp( [ 'Points remaining to be solved for: ' int2str(numPointsExtendedROI) ] )


%% Mapping IEP to Edge
if isMappingHarmonicNeighbourhoods
    %---
    % Determine nearest edgePoint to each internalExpansionPoint
    disp('Determining Expansion Capacity...')
    disp(['Using ' int2str( Options.numCPU ) ' processor(s)'])
    
    if matlabpool('size') == 0
        matlabpool( Options.numCPU ) ;
    end
    
    tic
    
    XYZedges = [ Options.voxelSize(1)*X1edge'; Options.voxelSize(2)*Y1edge'; Options.voxelSize(3)*Z1edge' ]; %
    
    %---
    % Restructure the way internalExpansionPoints are stored for the sake of the 'parfor' loop
    partXo   = cell(1, Options.numCPU);
    partYo   = partXo;
    partZo   = partXo;
    
    partIEPs      = floor( numIEPs/Options.numCPU ) ; % fraction of IEPs assigned to numCPU-1 processors
    lastPartIEPs  = numIEPs -( Options.numCPU - 1 )*partIEPs ; % largest portion, assigned to the last processor

            
    for cpu = 1 : Options.numCPU
        if cpu < Options.numCPU
            partXo{cpu} = Options.voxelSize(1)* Xo( (cpu-1)*partIEPs + 1 : cpu*partIEPs) ;
            partYo{cpu} = Options.voxelSize(2)* Yo( (cpu-1)*partIEPs + 1 : cpu*partIEPs) ;
            partZo{cpu} = Options.voxelSize(3)* Zo( (cpu-1)*partIEPs + 1 : cpu*partIEPs) ;
            
        elseif cpu == Options.numCPU
            partXo{cpu} = Options.voxelSize(1)* Xo( (cpu-1)*partIEPs + 1 : numIEPs) ;
            partYo{cpu} = Options.voxelSize(2)* Yo( (cpu-1)*partIEPs + 1 : numIEPs) ;
            partZo{cpu} = Options.voxelSize(3)* Zo( (cpu-1)*partIEPs + 1 : numIEPs) ;
        end
    end    
    
    partIEPtoEdgeDistance     = zeros( partIEPs, Options.numCPU ) ;
    lastPartIEPtoEdgeDistance = zeros( lastPartIEPs, Options.numCPU ) ;
    
    
        %---
        % Distance Calculations :
        % assign a portion of the IEPs to each available processor
        
        disp('Progress: ')
        parfor cpu = 1 : Options.numCPU
            
            XoCPU = partXo{cpu} ;
            YoCPU = partYo{cpu} ;
            ZoCPU = partZo{cpu} ;
            
            if cpu == Options.numCPU
                
                lastPartEdgeDist = zeros( lastPartIEPs, 1) ;
                
                for ko = 1 : lastPartIEPs
                    
                    disp( ['mapping... IEPtoEdge ' num2str(ko/lastPartIEPs, 2) ])

                    D                    = ( XYZedges - repmat( [ XoCPU(ko); YoCPU(ko); ZoCPU(ko) ], [ 1 numEdgePoints ] ) ) ;
                    lastPartEdgeDist(ko) = min( dot( D, D, 1 ) .^0.5 ) ;
                    
                end
                
                lastPartIEPtoEdgeDistance(:, cpu) = lastPartEdgeDist ;
                
            else
                
                partEdgeDist = zeros(partIEPs, 1) ;
                
                for ko = 1 : partIEPs
                    
                    D                = ( XYZedges - repmat( [ XoCPU(ko); YoCPU(ko); ZoCPU(ko) ], [ 1 numEdgePoints ] ) ) ;
                    partEdgeDist(ko) = min( dot( D, D, 1 ) .^0.5 ) ;
                                        
                end
                
                partIEPtoEdgeDistance(:, cpu) = partEdgeDist ;
            end
      
        end    
        
        lastPartIEPtoEdgeDistance = lastPartIEPtoEdgeDistance( :, Options.numCPU ) ;
        partIEPtoEdgeDistance     = partIEPtoEdgeDistance( :, 1 : Options.numCPU-1 ) ;
        
        IEPtoEdgeDistance         = [ partIEPtoEdgeDistance(:)' lastPartIEPtoEdgeDistance' ]' ;
        
        ComputationTime.mappingExpansionCapacity = toc ;
        clear tmp XYZedges partXo partYo partZo 
        
        if Options.isSavingInterimVar
        save( strcat(Options.name,'_IEPtoEdgeDistance'), 'IEPtoEdgeDistance') ;
        end
else

load( Options.IEPtoEdgeDistance ) ;

end



%% Mapping EP to IEP        
    if isMappingExpansionPoints
        %---
        % Determine nearest expansionPoints for voxels in the extendedROI
        disp('Determining best expansion points...')
        tic
        XYZo = [ Options.voxelSize(1)*Xo'; Options.voxelSize(2)*Yo'; Options.voxelSize(3)*Zo' ]; %
        
        if matlabpool('size') == 0
            matlabpool( Options.numCPU ) ;
        end
        
        disp(['Using ' int2str( Options.numCPU ) ' processor(s)'])
        
        partExtendedROI      = floor( numPointsExtendedROI/Options.numCPU ) ; % fraction of the extendedROI assigned to numCPU-1 processors
        lastPartExtendedROI  = numPointsExtendedROI -( Options.numCPU - 1 )*partExtendedROI ; % largest portion, assigned to the last processor
        
        %---
        % Restructure the way the extended ROI is stored for the sake of the 'parfor' loop
        partX1   = cell(1, Options.numCPU); 
        partY1   = partX1;
        partZ1   = partX1;
        
        for cpu = 1 : Options.numCPU
            if cpu < Options.numCPU
                partX1{cpu} = Options.voxelSize(1)* X1( (cpu-1)*partExtendedROI + 1 : cpu*partExtendedROI) ;
                partY1{cpu} = Options.voxelSize(2)* Y1( (cpu-1)*partExtendedROI + 1 : cpu*partExtendedROI) ;
                partZ1{cpu} = Options.voxelSize(3)* Z1( (cpu-1)*partExtendedROI + 1 : cpu*partExtendedROI) ;
                
            elseif cpu == Options.numCPU
                partX1{cpu} = Options.voxelSize(1)* X1( (cpu-1)*partExtendedROI + 1 : numPointsExtendedROI) ;
                partY1{cpu} = Options.voxelSize(2)* Y1( (cpu-1)*partExtendedROI + 1 : numPointsExtendedROI) ;
                partZ1{cpu} = Options.voxelSize(3)* Z1( (cpu-1)*partExtendedROI + 1 : numPointsExtendedROI) ;
            end
        end
        
        partIEPIndices            = zeros( partExtendedROI, Options.numCPU ) ;
        lastPartIEPIndices        = zeros( lastPartExtendedROI, Options.numCPU ) ;
        
        partIEPtoEPDistance       = zeros( partExtendedROI, Options.numCPU ) ;
        lastPartIEPtoEPDistance   = zeros( lastPartExtendedROI, Options.numCPU ) ;
        
        partExpansionCapacity     = zeros( partExtendedROI, Options.numCPU ) ;
        lastPartExpansionCapacity = zeros( lastPartExtendedROI, Options.numCPU ) ;
        
        %---
        % Distance Calculations :
        % assign a portion of the extended ROI to each available processor
        
        disp('Progress: ')
        parfor cpu = 1 : Options.numCPU
            
            X1CPU = partX1{cpu} ;
            Y1CPU = partY1{cpu} ;
            Z1CPU = partZ1{cpu} ;
            
            if cpu == Options.numCPU
                
                lastpartInd   = zeros( lastPartExtendedROI, 1) ;
                lastpartDist  = zeros( lastPartExtendedROI, 1) ;
                lastPartExCap = zeros( lastPartExtendedROI, 1) ;
                
                for k1 = 1 : lastPartExtendedROI
                    disp(['mapping... EPtoIEP ' num2str(k1/lastPartExtendedROI, 2)])

                    D        = repmat( [ X1CPU(k1); Y1CPU(k1); Z1CPU(k1) ], [ 1 numIEPs ] )-XYZo   ;
                    
                    [ IEPtoEPDist, idealIEP ]  = min( dot( D, D, 1 ) .^0.5 ) ;% Determines the nearest IEP to EP XYZ1(k)
                                        
                    lastpartDist(k1)  = IEPtoEPDist ;
                    lastpartInd(k1)   = sub2ind( gridDimensionVector, Xo(idealIEP), Yo(idealIEP), Zo(idealIEP) ) ;                
                    lastPartExCap(k1) = IEPtoEdgeDistance(idealIEP) ;
                    
                end
                
                lastPartIEPtoEPDistance(:, cpu)   = lastpartDist ;
                lastPartIEPIndices(:, cpu)        = lastpartInd ;              
                lastPartExpansionCapacity(:, cpu) = lastPartExCap ;
            else
                
                partInd   = zeros(partExtendedROI, 1) ;
                partDist  = zeros(partExtendedROI, 1) ;
                partExCap = zeros(partExtendedROI, 1) ;
                
                for k1 = 1 : partExtendedROI
                    D        = repmat([ X1CPU(k1); Y1CPU(k1); Z1CPU(k1) ], [ 1 numIEPs ] ) -XYZo ;
                    
                    [ IEPtoEPDist, idealIEP ]  = min( dot( D, D, 1) .^0.5 ) ;
                    partDist(k1)               = IEPtoEPDist ;
                    partInd(k1)                = sub2ind( gridDimensionVector, Xo(idealIEP), Yo(idealIEP), Zo(idealIEP) ) ;
                    partExCap(k1)              = IEPtoEdgeDistance(idealIEP) ;
                end
                
 
                partIEPtoEPDistance( :, cpu ) = partDist ;
                partIEPIndices( :, cpu )      = partInd ;
                partExpansionCapacity(:, cpu) = partExCap ;
            end
                        
        end           
        
        matlabpool close
        
        lastPartIEPIndices          = lastPartIEPIndices( :, Options.numCPU ) ;
        partIEPIndices              = partIEPIndices( :, 1 : Options.numCPU-1 ) ;
        
        lastPartIEPtoEPDistance     = lastPartIEPtoEPDistance( :, Options.numCPU ) ;
        partIEPtoEPDistance         = partIEPtoEPDistance( :, 1 : Options.numCPU-1 ) ;
        
        lastPartExpansionCapacity   = lastPartExpansionCapacity( :, Options.numCPU ) ;
        partExpansionCapacity       = partExpansionCapacity( :, 1 : Options.numCPU-1 ) ;
        
        associatedIEPs              = [ partIEPIndices(:)' lastPartIEPIndices' ]' ; % Expansion points (linear indices paired w/points in extendedROI)
        associatedIEPtoEPDistance   = [ partIEPtoEPDistance(:)' lastPartIEPtoEPDistance' ]' ; % SO WHAT THE NAME IS REDUNDANT W/E
        associatedIEPtoEdgeDistance = [ partExpansionCapacity(:)' lastPartExpansionCapacity' ]' ;
        
        ComputationTime.mappingIEPs = toc
        
        if Options.isSavingInterimVar
        save( strcat(Options.name, '_IEPtoEPassociations'), 'associatedIEPs', 'associatedIEPtoEPDistance', 'associatedIEPtoEdgeDistance') ;   
        end
        
    else
        
        load( Options.IEPtoEPassociations ) ;
        
    end
    

    

%% gradient calculations
if isCalculatingGradients
    disp('Calculating gradients')
    tic    
        
    GradientTerms = sparsegradient( reducedBackgroundField, expansionSurface, Options.expansionOrder );
    
    ComputationTime.gradients = toc 
        
    if Options.isSavingInterimVar 
        save( strcat(Options.name,'_GradientTerms'), 'GradientTerms' )
    end
    
else
    load( Options.GradientTerms ) ;
end
     
    
    %% Taylor expansion
    
    % locations of the IEPs for each EP...
    [Xiep, Yiep, Ziep] = ind2sub( gridDimensionVector, squeeze(associatedIEPs) ) ;
    
    % not a vector but a set of column vectors...
    distanceVector = [  (X1 - Xiep)';   (Y1 - Yiep)' ;  (Z1 - Ziep)' ] ; 
    
    if isPerformingTaylorExpansion
        %---
        % Estimate extendedBackgroundField...
        disp('Performing Taylor expansion...')
        tic
        
        % Housekeeping/rearrangement...
        GradientTerms.directions = permute( GradientTerms.directions, [ 1 3 2 ] ) ;
        
        % Physically meaningless but simplifies computation...
        complementGDirections    = GradientTerms.directions == 0 ;
        
        % couldn't think of a better name...
        factorialFactor          = factorial( sum( sum( GradientTerms.directions, 1 ), 3 ) ) ;
        
        % series terms
        taylorTerms              = zeros( [ numPointsExtendedROI (Options.expansionOrder + 1) ] ) ;
        
        for k1 = 1 : numPointsExtendedROI
            
            disp( ['extrapolating...' num2str(k1/numPointsExtendedROI, 2)] ) ;
            
            ko = GradientTerms.indexKey( associatedIEPs(k1) ) ; % IEP ko corresponding to external extension point k1
            
            % the "dx*dx*dy*dy*dz"-type factor:
            tmp =          distanceVector(1, k1)*GradientTerms.directions(:, :, 1) + complementGDirections(:, :, 1) ;
            tmp = tmp .* ( distanceVector(2, k1)*GradientTerms.directions(:, :, 2) + complementGDirections(:, :, 2) ) ;
            tmp = tmp .* ( distanceVector(3, k1)*GradientTerms.directions(:, :, 3) + complementGDirections(:, :, 3) ) ;
            tmp = prod(tmp, 1) ./ factorialFactor ;
            
            % 0th order Taylor:
            taylorTerms(k1, 1) = GradientTerms.coefficients(ko, 1) ;
            
            % Higher orders
            for order = 2 : Options.expansionOrder + 1
                taylorTerms(k1, order) = sum( GradientTerms.coefficients(ko, gradInd(order, 1):gradInd(order, 2)) .* tmp(1, gradInd(order, 1):gradInd(order, 2))) ;
            end
            
        end
        
        if Options.isSavingInterimVar
            save( strcat(Options.name, '_taylorTerms'), 'taylorTerms' ) ;
        end
    
    else
        load( Options.taylorTerms )
    end
    
    % combine individual series terms
    taylorSeries = cumsum( taylorTerms, 2) ;
    ComputationTime.taylor = toc
    
    % reshape
    extendedBackgroundField = zeros( [ gridDimensionVector Options.expansionOrder + 1] );
    
    for order = 1 : Options.expansionOrder + 1
        tmp                                   = zeros( gridDimensionVector ) ;
        tmp( extensionPoints )                = taylorSeries(:, order);
        extendedBackgroundField(:,:,:, order) = tmp ;
        
        tmp = [];
    end
    
        clear taylorTerms taylorSeries


%% Reformat...
        
        expansionSurface = reshape(expansionSurface, gridDimensionVector) ;
        
        tmp                            = expansionSurface ;
        tmp( internalExpansionPoints ) = IEPtoEdgeDistance ;
        IEPtoEdgeDistance              = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = associatedIEPtoEdgeDistance ;
        associatedIEPtoEdgeDistance    = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = associatedIEPs ;
        associatedIEPs                 = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = associatedIEPtoEPDistance ;
        associatedIEPtoEPDistance      = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = distanceVector(1,:) ;
        IEPtoEPDistanceXYZ             = tmp ;
        tmp( extensionPoints )         = distanceVector(2,:) ;
        IEPtoEPDistanceXYZ(:,:,:,2)    = tmp ;
        tmp( extensionPoints )         = distanceVector(3,:) ;
        IEPtoEPDistanceXYZ(:,:,:,3)    = tmp ;
        
        % if the EP is further from the IEP than the nearest egde, it is
        % outside the harmonic neighbourhood and it will have to be solved
        % next iteration w/a different IEP
        %
        % possible issue w/ coding:
        % assumption is that the distance vector from EP to nearest IEP won't cross
        % ROI boundaries. This should be OK for the brain but may pose a problem for other geometries.
        
        expansionCapacity = abs( ceil(associatedIEPtoEPDistance) ) < abs( floor(associatedIEPtoEdgeDistance) ) ;
        
    
    
    %% Assess convergence
    diffPenultimateOrder = abs( extendedBackgroundField(:,:,:, Options.expansionOrder - 1) - extendedBackgroundField(:,:,:, Options.expansionOrder ) );
    diffHighestOrder     = abs( extendedBackgroundField(:,:,:, Options.expansionOrder) - extendedBackgroundField(:,:,:, Options.expansionOrder + 1) );
    
    convergenceSurface   =  diffPenultimateOrder > diffHighestOrder;
    
    clear diffHighestOrder diffPenultimateOrder
    
    
    %% Output
    
    EdgeOut.extendedBackgroundField = extendedBackgroundField ;
    extendedBackgroundField         = extendedBackgroundField(:,:,:, Options.expansionOrder + 1) ;
    
    EdgeOut.extendedROI             = extendedROI ; 
    clear extendedROI
        
    EdgeOut.numPointsExtendedROI    = numPointsExtendedROI ;
    clear numPointsExtendedROI
       
    if currentIteration == 1
        EdgeOut.reducedBackgroundField  = inputField ;
    else
        EdgeOut.reducedBackgroundField  = reducedBackgroundField ;
    end
    
    EdgeOut.IEPtoEdgeDistance         = IEPtoEdgeDistance ;
    clear IEPtoEdgeDistance
    
    EdgeOut.associatedIEPtoEPDistance = associatedIEPtoEPDistance ; % Euclidean distance
    clear associatedIEPtoEPDistance 

    EdgeOut.IEPtoEPDistanceXYZ        = IEPtoEPDistanceXYZ ; % 4D Array - distance in voxels
    clear IEPtoEPDistanceXYZ
        
    EdgeOut.associatedIEPs            = associatedIEPs ;
    clear associatedIEPs 
    
    EdgeOut.ROI                     = ROI ;
    EdgeOut.reducedROI              = reducedROI ;
    EdgeOut.convergenceSurface      = convergenceSurface ;
    EdgeOut.expansionCapacity       = expansionCapacity ;

    
    ComputationTime.total           = toc ;
    EdgeOut.ComputationTime         = ComputationTime ;
    
    EdgeOut.Options                 = Options ;
    
    if Options.isWritingToDisk
    disp('Saving to disk...')
    save( Options.name, 'EdgeOut', '-v7.3')
    end
    
    if currentIteration < Options.numIterations
        clear EdgeOut
    end
        
    %% Update
    isMappingExpansionPoints        = true ;
    isMappingHarmonicNeighbourhoods = true ;
    isCalculatingGradients          = true ;
    isPerformingTaylorExpansion     = true ;
    
    if Options.isUsingConvergenceCriterion
        reducedROI         = reducedROI + convergenceSurface .* expansionCapacity ;
    else
        reducedROI         = reducedROI + expansionCapacity ;
    end
    
    reducedBackgroundField = reducedBackgroundField + reducedROI .* extendedBackgroundField ;
    
    previousIEPs           = previousIEPs + expansionSurface ;
        
end


end