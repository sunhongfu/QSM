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
%           iteration using different IPs. 
%               default: false
%              
%       isDisplayingProgress
%           incremental progress printed to screen
%           default: false
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
%       IPtoEdgeDistance
%       specifies path of MATLAB matrix containing the variables
%           -IPtoEdgeDistance : determines 'harmonic neighbourhood' about
%           each IP
%               
%       IPtoEPAssociations
%           specifies path of MATLAB matrix containing the variables
%           -associatedIPs : linear indices of IPs for each EP 
%           
%           -associatedIPtoEPDistance : Euclidean distance b/tw them
%           
%           -associatedIPtoEdgeDistance : Each EP contains the
%           IPtoEdgeDistance of its corresponding IP.               
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
%       IP = Internal Point : point at which extrapolation is based   
%       EP  = Extension Point : looking to determine Bkgr. Field here 
%
%
% R Topfer 2012     topfer@ualberta.ca
%

% TODO
% Easy final step:
% Examine points outside 'harmonic neighbourhood'
% Rather than perform 2nd iteration, simply assign them the bkgr field of
% the nearest valid point (i.e. 0th order expansion)
%
% Even/Odd dimensions compatibility


%% constants

DEFAULT_NAME                            = 'sharpEdgesResults' ;
DEFAULT_EXPANSIONORDER                  = 2 ;
DEFAULT_NUMITERATIONS                   = 1 ;
DEFAULT_NUMCPU                          = parcluster ;
DEFAULT_NUMCPU                          = DEFAULT_NUMCPU.NumWorkers ;
DEFAULT_VOXELSIZE                       = [1 1 1] ;

DEFAULT_ISUSINGCONVERGENCECRITERION     = false ;

DEFAULT_ISMAPPINGEXPANSIONPOINTS        = true ;
DEFAULT_ISMAPPINGHARMONICNEIGHBOURHOODS = true ;
DEFAULT_ISCALCULATINGGRADIENTS          = true ;
DEFAULT_ISPERFORMINGTAYLOREXPANSION     = true ;

DEFAULT_ISWRITINGTODISK                 = true ;
DEFAULT_ISDISPLAYINGPROGRESS            = false ;


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
    disp('Options.voxelSize not assigned. Assuming isotropic [1 1 1].');
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end

if ~myisfield( Options, 'aspectRatio' ) || isempty(Options.aspectRatio)
    Options.aspectRatio = round( Options.voxelSize / min( Options.voxelSize ) ) ;
end

if ~myisfield( Options, 'isUsingConvergenceCriterion') || isempty(Options.isUsingConvergenceCriterion)
    Options.isUsingConvergenceCriterion = DEFAULT_ISUSINGCONVERGENCECRITERION ;
end

if  ~myisfield( Options, 'IPtoEPAssociations' ) || isempty(Options.IPtoEPAssociations)
    isMappingExpansionPoints = DEFAULT_ISMAPPINGEXPANSIONPOINTS ;
else
    isMappingExpansionPoints = false ;
end

if ~myisfield( Options, 'IPtoEdgeDistance') || isempty(Options.IPtoEdgeDistance)
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

if ~myisfield( Options, 'isDisplayingProgress') || isempty(Options.isDisplayingProgress)
    Options.isDisplayingProgress = DEFAULT_ISDISPLAYINGPROGRESS ;
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
previousIPs            = zeros( gridDimensionVector ) ;

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
    tmpROIedge       = tmpROI - shaver(tmpROI, 1 ) ; % avoid using the same IPs one iteration to the next...                      
    tmpROI           = tmpROI - tmpROIedge .* previousIPs ; % concern is exclusively w/previous IPs on tmpROIedge : w/out this, in the next line, we'd introduce more pt.s into the pool of potential IPs that would not be chosen/used but their inclusion would slow down mapping steps  
    tmpROI           = tmpROI - shaver(tmpROI, 1 ) ; % reduced to ~ lamina
    expansionSurface = tmpROI(:) ;
    
    clear tmpROI tmpROIedge
    
    internalExpansionPoints = find( expansionSurface ) ;
    numIPs                 = length( internalExpansionPoints ) ; % num of pt.s to choose from for determining best expansion point
    [Xo, Yo, Zo]            = ind2sub( gridDimensionVector, internalExpansionPoints ) ;
    
    if Options.isDisplayingProgress
        
        disp( [ 'Solving for extended background field: Iteration ' int2str(currentIteration) ])
        disp( [ 'Points remaining to be solved for: ' int2str(numPointsExtendedROI) ] )
        disp(['Using ' int2str( Options.numCPU ) ' processor(s)'])
    end
    


    
%% Mapping IP to Edge
if isMappingHarmonicNeighbourhoods
        %---
        % Determine nearest edgePoint to each internalPoint
        disp('Determining Expansion Capacity...')
    
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
    
    partIPs      = floor( numIPs/Options.numCPU ) ; % fraction of IPs assigned to numCPU-1 processors
    lastPartIPs  = numIPs -( Options.numCPU - 1 )*partIPs ; % largest portion, assigned to the last processor

            
    for cpu = 1 : Options.numCPU
        if cpu < Options.numCPU
            partXo{cpu} = Options.voxelSize(1)* Xo( (cpu-1)*partIPs + 1 : cpu*partIPs) ;
            partYo{cpu} = Options.voxelSize(2)* Yo( (cpu-1)*partIPs + 1 : cpu*partIPs) ;
            partZo{cpu} = Options.voxelSize(3)* Zo( (cpu-1)*partIPs + 1 : cpu*partIPs) ;
            
        elseif cpu == Options.numCPU
            partXo{cpu} = Options.voxelSize(1)* Xo( (cpu-1)*partIPs + 1 : numIPs) ;
            partYo{cpu} = Options.voxelSize(2)* Yo( (cpu-1)*partIPs + 1 : numIPs) ;
            partZo{cpu} = Options.voxelSize(3)* Zo( (cpu-1)*partIPs + 1 : numIPs) ;
        end
    end    
    
    partIPtoEdgeDistance     = zeros( partIPs, Options.numCPU ) ;
    lastPartIPtoEdgeDistance = zeros( lastPartIPs, Options.numCPU ) ;
    
    
        %---
        % Distance Calculations :
        % assign a portion of the IPs to each available processor
        
        parfor cpu = 1 : Options.numCPU
            
            XoCPU = partXo{cpu} ;
            YoCPU = partYo{cpu} ;
            ZoCPU = partZo{cpu} ;
            
            if cpu == Options.numCPU
                
                lastPartEdgeDist = zeros( lastPartIPs, 1) ;
                if Options.isDisplayingProgress
                    for ko = 1 : lastPartIPs
                        
                        disp( ['mapping... IPtoEdge ' num2str(ko/lastPartIPs, 2) ])
                        
                        D                    = ( XYZedges - repmat( [ XoCPU(ko); YoCPU(ko); ZoCPU(ko) ], [ 1 numEdgePoints ] ) ) ;
                        lastPartEdgeDist(ko) = min( dot( D, D, 1 ) .^0.5 ) ;
                        
                    end
                    
                else
                    for ko = 1 : lastPartIPs
                        
                        D                    = ( XYZedges - repmat( [ XoCPU(ko); YoCPU(ko); ZoCPU(ko) ], [ 1 numEdgePoints ] ) ) ;
                        lastPartEdgeDist(ko) = min( dot( D, D, 1 ) .^0.5 ) ;
                        
                    end
                    lastPartIPtoEdgeDistance(:, cpu) = lastPartEdgeDist ;
                end
                
            else
                
                partEdgeDist = zeros(partIPs, 1) ;
                
                for ko = 1 : partIPs
                    
                    D                = ( XYZedges - repmat( [ XoCPU(ko); YoCPU(ko); ZoCPU(ko) ], [ 1 numEdgePoints ] ) ) ;
                    partEdgeDist(ko) = min( dot( D, D, 1 ) .^0.5 ) ;
                                        
                end
                
                partIPtoEdgeDistance(:, cpu) = partEdgeDist ;
            end
      
        end    
        
        lastPartIPtoEdgeDistance = lastPartIPtoEdgeDistance( :, Options.numCPU ) ;
        partIPtoEdgeDistance     = partIPtoEdgeDistance( :, 1 : Options.numCPU-1 ) ;
        
        IPtoEdgeDistance         = [ partIPtoEdgeDistance(:)' lastPartIPtoEdgeDistance' ]' ;
        
        ComputationTime.mappingExpansionCapacity = toc ;
        clear tmp XYZedges partXo partYo partZo 
        
        if Options.isSavingInterimVar
        save( strcat(Options.name,'_IPtoEdgeDistance'), 'IPtoEdgeDistance') ;
        end
else

load( Options.IPtoEdgeDistance ) ;

end



%% Mapping EP to IP        
    if isMappingExpansionPoints
        %---
        % Determine nearest expansionPoints for voxels in the extendedROI
        
        tic
        
        if matlabpool('size') == 0
            matlabpool( Options.numCPU ) ;
        end
        
        disp('Determining best expansion points...')
        
        XYZo = [ Options.voxelSize(1)*Xo'; Options.voxelSize(2)*Yo'; Options.voxelSize(3)*Zo' ]; %
        
        
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
        
        partIPIndices            = zeros( partExtendedROI, Options.numCPU ) ;
        lastPartIPIndices        = zeros( lastPartExtendedROI, Options.numCPU ) ;
        
        partIPtoEPDistance       = zeros( partExtendedROI, Options.numCPU ) ;
        lastPartIPtoEPDistance   = zeros( lastPartExtendedROI, Options.numCPU ) ;
        
        partExpansionCapacity     = zeros( partExtendedROI, Options.numCPU ) ;
        lastPartExpansionCapacity = zeros( lastPartExtendedROI, Options.numCPU ) ;
        
        %---
        % Distance Calculations :
        % assign a portion of the extended ROI to each available processor

        parfor cpu = 1 : Options.numCPU
            
            X1CPU = partX1{cpu} ;
            Y1CPU = partY1{cpu} ;
            Z1CPU = partZ1{cpu} ;
            
            if cpu == Options.numCPU
                
                lastpartInd   = zeros( lastPartExtendedROI, 1) ;
                lastpartDist  = zeros( lastPartExtendedROI, 1) ;
                lastPartExCap = zeros( lastPartExtendedROI, 1) ;
                
                if Options.isDisplayingProgress
                    for k1 = 1 : lastPartExtendedROI
                        
                        disp(['mapping... EPtoIP ' num2str(k1/lastPartExtendedROI, 2)])
                        
                        D        = repmat( [ X1CPU(k1); Y1CPU(k1); Z1CPU(k1) ], [ 1 numIPs ] )-XYZo   ;
                        
                        [ IPtoEPDist, idealIP ]  = min( dot( D, D, 1 ) .^0.5 ) ;% Determines the nearest IP to EP XYZ1(k)
                        
                        lastpartDist(k1)  = IPtoEPDist ;
                        lastpartInd(k1)   = sub2ind( gridDimensionVector, Xo(idealIP), Yo(idealIP), Zo(idealIP) ) ;
                        lastPartExCap(k1) = IPtoEdgeDistance(idealIP) ;
                        
                    end
                    
                else
                    
                    for k1 = 1 : lastPartExtendedROI

                        D        = repmat( [ X1CPU(k1); Y1CPU(k1); Z1CPU(k1) ], [ 1 numIPs ] )-XYZo   ;
                        
                        [ IPtoEPDist, idealIP ]  = min( dot( D, D, 1 ) .^0.5 ) ;% Determines the nearest IP to EP XYZ1(k)
                        
                        lastpartDist(k1)  = IPtoEPDist ;
                        lastpartInd(k1)   = sub2ind( gridDimensionVector, Xo(idealIP), Yo(idealIP), Zo(idealIP) ) ;
                        lastPartExCap(k1) = IPtoEdgeDistance(idealIP) ;
                        
                    end
                end
                lastPartIPtoEPDistance(:, cpu)   = lastpartDist ;
                lastPartIPIndices(:, cpu)        = lastpartInd ;              
                lastPartExpansionCapacity(:, cpu) = lastPartExCap ;
            else
                
                partInd   = zeros(partExtendedROI, 1) ;
                partDist  = zeros(partExtendedROI, 1) ;
                partExCap = zeros(partExtendedROI, 1) ;
                
                for k1 = 1 : partExtendedROI
                    D        = repmat([ X1CPU(k1); Y1CPU(k1); Z1CPU(k1) ], [ 1 numIPs ] ) -XYZo ;
                    
                    [ IPtoEPDist, idealIP ]  = min( dot( D, D, 1) .^0.5 ) ;
                    partDist(k1)               = IPtoEPDist ;
                    partInd(k1)                = sub2ind( gridDimensionVector, Xo(idealIP), Yo(idealIP), Zo(idealIP) ) ;
                    partExCap(k1)              = IPtoEdgeDistance(idealIP) ;
                end
                
 
                partIPtoEPDistance( :, cpu ) = partDist ;
                partIPIndices( :, cpu )      = partInd ;
                partExpansionCapacity(:, cpu) = partExCap ;
            end
                        
        end           
        
        matlabpool close
        
        lastPartIPIndices          = lastPartIPIndices( :, Options.numCPU ) ;
        partIPIndices              = partIPIndices( :, 1 : Options.numCPU-1 ) ;
        
        lastPartIPtoEPDistance     = lastPartIPtoEPDistance( :, Options.numCPU ) ;
        partIPtoEPDistance         = partIPtoEPDistance( :, 1 : Options.numCPU-1 ) ;
        
        lastPartExpansionCapacity   = lastPartExpansionCapacity( :, Options.numCPU ) ;
        partExpansionCapacity       = partExpansionCapacity( :, 1 : Options.numCPU-1 ) ;
        
        associatedIPs              = [ partIPIndices(:)' lastPartIPIndices' ]' ; % Expansion points (linear indices paired w/points in extendedROI)
        associatedIPtoEPDistance   = [ partIPtoEPDistance(:)' lastPartIPtoEPDistance' ]' ; % SO WHAT THE NAME IS REDUNDANT W/E
        associatedIPtoEdgeDistance = [ partExpansionCapacity(:)' lastPartExpansionCapacity' ]' ;
        
        ComputationTime.mappingIPs = toc ;
        
        if Options.isSavingInterimVar
        save( strcat(Options.name, '_IPtoEPAssociations'), 'associatedIPs', 'associatedIPtoEPDistance', 'associatedIPtoEdgeDistance') ;   
        end
        
    else
        
        load( Options.IPtoEPAssociations ) ;
        
    end
    

    

%% gradient calculations
if isCalculatingGradients
    disp('Calculating gradients...')
    tic    
        
    GradientTerms = sparsegradient( reducedBackgroundField, expansionSurface, Options.expansionOrder );
    
    ComputationTime.gradients = toc ;
        
    if Options.isSavingInterimVar 
        save( strcat(Options.name,'_GradientTerms'), 'GradientTerms' )
    end
    
else
    load( Options.GradientTerms ) ;
end
     
    
    %% Taylor expansion
    
    % locations of the IPs for each EP...
    [Xiep, Yiep, Ziep] = ind2sub( gridDimensionVector, squeeze(associatedIPs) ) ;
    
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
        
                        if Options.isDisplayingProgress
        
        
                            for k1 = 1 : numPointsExtendedROI
                                
                                disp( ['extrapolating...' num2str(k1/numPointsExtendedROI, 2)] ) ;
                                
                                ko = GradientTerms.indexKey( associatedIPs(k1) ) ; % IP ko corresponding to external extension point k1
                                
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
                        else
                            for k1 = 1 : numPointsExtendedROI
                                
                                ko = GradientTerms.indexKey( associatedIPs(k1) ) ; % IP ko corresponding to external extension point k1
                                
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
                        end
        
        if Options.isSavingInterimVar
            save( strcat(Options.name, '_taylorTerms'), 'taylorTerms' ) ;
        end
    
    else
        load( Options.taylorTerms )
    end
    
    % combine individual series terms
    taylorSeries = cumsum( taylorTerms, 2) ;
    ComputationTime.taylor = toc ;
    
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
        tmp( internalExpansionPoints ) = IPtoEdgeDistance ;
        IPtoEdgeDistance              = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = associatedIPtoEdgeDistance ;
        associatedIPtoEdgeDistance    = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = associatedIPs ;
        associatedIPs                 = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = associatedIPtoEPDistance ;
        associatedIPtoEPDistance      = tmp ;
        
        tmp                            = extendedROI ;
        tmp( extensionPoints )         = distanceVector(1,:) ;
        IPtoEPDistanceXYZ             = tmp ;
        tmp( extensionPoints )         = distanceVector(2,:) ;
        IPtoEPDistanceXYZ(:,:,:,2)    = tmp ;
        tmp( extensionPoints )         = distanceVector(3,:) ;
        IPtoEPDistanceXYZ(:,:,:,3)    = tmp ;
        
        % if the EP is further from the IP than the nearest egde, it is
        % outside the harmonic neighbourhood and it will have to be solved
        % next iteration w/a different IP
        %
        % possible issue w/ coding:
        % assumption is that the distance vector from EP to nearest IP won't cross
        % ROI boundaries. This should be OK for the brain but may pose a problem for other geometries.
        
        expansionCapacity = abs( ceil(associatedIPtoEPDistance) ) < abs( floor(associatedIPtoEdgeDistance) ) ;
        
    
    
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
    
    EdgeOut.IPtoEdgeDistance         = IPtoEdgeDistance ;
    clear IPtoEdgeDistance
    
    EdgeOut.associatedIPtoEPDistance = associatedIPtoEPDistance ; % Euclidean distance
    clear associatedIPtoEPDistance 

    EdgeOut.IPtoEPDistanceXYZ        = IPtoEPDistanceXYZ ; % 4D Array - distance in voxels
    clear IPtoEPDistanceXYZ
        
    EdgeOut.associatedIPs            = associatedIPs ;
    clear associatedIPs 
    
    EdgeOut.ROI                     = ROI ;
    EdgeOut.reducedROI              = reducedROI ;
    EdgeOut.convergenceSurface      = convergenceSurface ;
    EdgeOut.expansionCapacity       = expansionCapacity ;

    
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
    
    previousIPs           = previousIPs + expansionSurface ;
        
end


end