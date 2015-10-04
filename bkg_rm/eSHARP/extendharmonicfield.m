function[ harmonicFieldExtended, A, M ] = extendharmonicfield ...
	( harmonicFieldReduced, mask, maskReduced, Parameters )
%EXTENDHARMONICFIELD recover/extrapolate field 
%
%   Uses a finite-difference Taylor expansion to extrapolate a harmonic 
%   field map (e.g. "background" field over a reduced inner portion of the 
%   brain) to an extended region (e.g. full brain up to the cortex)
%
% =========================================================================
%
%               USAGE
%       .......................
%
%       Syntax
%   
%       [xEP,A,M] = EXTENDHARMONICFIELD( xIP, mask, maskIp, Parameters )    
%       
%       .......................
%
%       Description
%
%       EXTENDHARMONICFIELD takes a harmonic field (xIP) defined over some
%       reduced region of "internal points" (IP), themselves defined by 3D
%       binary array maskIp, and extrapolates the field out to a set of
%       "external points" (EP): viz, the voxels contained within the binary
%       array 'mask' which do not overlap with the IP of maskIp.
%       
%       The extrapolated field xEP is the product of M' * A * xIP(:) - reshaped
%       to the original dimensions of xIP - where M is a trunctation (masking)
%       matrix operator for the EP (such that M*x(:) 'picks out' the EP from
%       vector x), and A is a matrix operator consisting of Taylor-expansion
%       based coefficients. (In most cases, the user probably only cares about
%       xEP.)
%
%       Parameters must possess the following fields
%
%		.voxelSize
%			(e.g. .Parameters.voxelSize = [1 1 1] mm, isotropic)
%
%       .radius
%           determins size of spherical region over which extrapolation takes 
%           place.
%           ***NB*** if this is not large enough to incorporate all of the EP
%           the extrapolation will FAIL on account of the truncation operator
%           M being incorrect. (to be fixed.)
%
%       .......................
%
%
%       The following Parameters.fields are optional 
%
%             
%           .expansionOrder
%               default: 2
%
%           .isDisplayingProgress
%               default: true
%
%           .tmpSaveFldr
%               default: ./
%               (user must have write permission to folder)
%
%
%
% ABBRV:
%       IP = Internal Point 
%			 e.g. point at which extrapolation is based   
%
%       EP = Extension/Edge Point 
%			 e.g. point to which the field will be extrapolated
%
% =========================================================================
%       
%             REFERENCES
%       .......................
%
% Topfer et al. SHARP edges: Recovering cortical phase contrast through
% harmonic extension. Magn Reson Med. 2015;73(2):851-6. doi 10.1002/mrm.25148
%
% =========================================================================
%       
%               LICENSE
%       .......................
%   
%   Unlicensed, confidential, top secret. 
%                 ;)
%			  
%       2014 - topfer@ualberta.ca
% =========================================================================
%               
%               INSTALL               
%       .......................
% 
% *NB:
%   Extendharmonicfield() calls a number of miscellaneous Matlab functions,
% the scripts for which will of course need to be discoverable by Matlab 
% (i.e. in its search path). In addition, a call is made to a command line 
% tool: 
%       mapip2ep
%   
%   User must first BUILD the mapip2ep executable from mapipto2ep.cpp
%
%   User is free to put the resulting mapip2ep binary anywhere but must
%   inform Matlab of its location by changing the following local variable
%   within extendharmonicfield.m: 
%
%       pathToBinaries = '../bin/'
%
% =========================================================================
% =========================================================================

pathToBinaries = '~/Projects/General/bin/' ;

% =========================================================================
% =========================================================================


% -------------------------------------------------------------------------
% Add mapip2ep to the path
setenv('PATH', [getenv('PATH') ':' pathToBinaries]) ;

% -------
% Check mapip2ep is within the path
[returnedStatus,~] = system('mapip2ep ') ;
if returnedStatus ~= 1
    disp('Error: see *NB in function HELP:')
    help(mfilename);
    return;
end

% -------------------------------------------------------------------------
% Examine input arguments 
if nargin < 4 || isempty(Parameters) 
    disp('Error: function requires 4 arguments')
    help(mfilename);
    return;
end

DEFAULT_EXPANSIONORDER          = 2 ;
DEFAULT_ISDISPLAYINGPROGRESS    = true ; 
DEFAULT_TMPSAVEFLDR             = './' ;


if  ~myisfield( Parameters, 'voxelSize' ) || (nnz(Parameters.voxelSize) ~=3)
    error('Parameters.voxelSize must be a 3 component vector. See HELP') ;
end

if  ~myisfield( Parameters, 'radius' ) || isempty(Parameters.radius) 
    error('Parameters.radius must be defined (a single or 3 component vector). See HELP') ;
elseif length(Parameters.radius) == 1 
    Parameters.radius = Parameters.radius*[1 1 1] ;
end

if  ~myisfield( Parameters, 'expansionOrder' ) || isempty(Parameters.expansionOrder)
    Parameters.expansionOrder = DEFAULT_EXPANSIONORDER ;
end

if  ~myisfield( Parameters, 'isDisplayingProgress' ) || isempty(Parameters.isDisplayingProgress)
    Parameters.isDisplayingProgress = DEFAULT_ISDISPLAYINGPROGRESS ;
end

if  ~myisfield( Parameters, 'tmpSaveFldr' ) || isempty(Parameters.tmpSaveFldr)
    Parameters.tmpSaveFldr = DEFAULT_TMPSAVEFLDR ;
end



% How 'far back' (# voxels) from the edge of the maskReduced the IP are
Parameters.offset = Parameters.expansionOrder ; 


% =========================================================================

% -------
% Create extension operator

maskEp = logical( mask - maskReduced );

% outerEdgeEP = mask - shaver( mask, 1 ) ;
% outerEdgeIP = maskReduced - shaver( maskReduced, 1 ) ;
%
%     
% % -------------------------------------------------------------------------
%
% % -------------------------------------------------------------------------
% % To determine the radius of the the harmonic neighbourhood overwhich the
% % expansion takes place, the following will determine the max distance 
% % between the outermost EP and the outermost IP.
% % The radius will then be set to 
% %                               ceil( (max.dist + 1 )/2) 
% %
% % This ensures all EP will be covered, but it isn't necessarily faithful 
% % to the underlying physics (i.e. it implicitly assumes that the spherical
% % region centered at IP and of radius sqrt((EP-IP)^2) does NOT include any
% % sources of background field (e.g. an air-tissue interface)). 
%
%
% [xOuterEdgeEP, yOuterEdgeEP, zOuterEdgeEP] = ...
%     ind2sub( size(outerEdgeEP), find( outerEdgeEP ) ) ;
%
% xOuterEdgeEP = Parameters.voxelSize(1)*xOuterEdgeEP;
% yOuterEdgeEP = Parameters.voxelSize(2)*yOuterEdgeEP;
% zOuterEdgeEP = Parameters.voxelSize(3)*zOuterEdgeEP;
%
% nOuterEdgeEP = length( xOuterEdgeEP ) ;
%
%
% [xOuterEdgeIP, yOuterEdgeIP, zOuterEdgeIP] = ...
%     ind2sub( size(outerEdgeIP), find( outerEdgeIP ) ) ;
%
% xyzOuterEdgeIP = [ Parameters.voxelSize(1)*xOuterEdgeIP'; 
%                    Parameters.voxelSize(2)*yOuterEdgeIP'; 
%                    Parameters.voxelSize(3)*zOuterEdgeIP' ]; 
% disp( ' size ' )
% size(xyzOuterEdgeIP)
%
% nOuterEdgeIP = length( xOuterEdgeIP ) ;
%
%
% maxMinIpToEpDistance = 0 ;
% minDistancesIPtoEP = zeros( nOuterEdgeEP, 1 ) ;
%
% if Parameters.isDisplayingProgress 
%     for iOuterEdgeEP = 1 : nOuterEdgeEP 
%         disp(num2str(100*iOuterEdgeEP/nOuterEdgeEP) ) ; 
%         distancesIPtoEP = xyzOuterEdgeIP - repmat( ...
%                     [ xOuterEdgeEP(iOuterEdgeEP); 
%                       yOuterEdgeEP(iOuterEdgeEP); 
%                       zOuterEdgeEP(iOuterEdgeEP); ], [ 1 nOuterEdgeIP ] ) ; 
%
%         minDistancesIPtoEP( iOuterEdgeEP ) = ...
%             min( dot( distancesIPtoEP, distancesIPtoEP, 1) .^ 0.5 ) ;
%         
%     end
% else
%     for iOuterEdgeEP = 1 : nOuterEdgeEP 
%         distancesIPtoEP = xyzOuterEdgeIP - repmat( ...
%                     [ xOuterEdgeEP(iOuterEdgeEP); 
%                       yOuterEdgeEP(iOuterEdgeEP); 
%                       zOuterEdgeEP(iOuterEdgeEP); ], [ 1 nOuterEdgeIP ] ) ; 
%
%         minDistancesIPtoEP( iOuterEdgeEP ) = ...
%             min( dot( distancesIPtoEP, distancesIPtoEP, 1) .^ 0.5 ) ;
%         
%     end
% end
%
% Parameters.radius = ceil( ceil(max( minDistancesIPtoEP )) ./Parameters.voxelSize ) ;


maskIp = logical( shaver( maskReduced, Parameters.offset ) ...
                - shaver( maskReduced, Parameters.radius ) );

[ A, M ] = createextensionoperator( maskIp, maskEp, Parameters ) ;

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
harmonicFieldExtended = reshape( M'*A*harmonicFieldReduced(:), ...
                                size(harmonicFieldReduced) ) ;


end



function[ A, M ] = createextensionoperator( maskIp, maskEp, Parameters )
%CREATEEXTENSIONOPERATOR
%
%	Returns a matrix operator to extrapolate a harmonic field
%
% =========================================================================
%
%       Syntax
%   
%       [A,M] = CREATEEXTENSIONOPERATOR( maskIp, maskEp, Parameters )
%       
%       .......................
%
%       Description
%
%       CREATEEXTENSIONOPERATOR takes a 3D binary array of "internal points"
%       (maskIp) and outputs an extension operator (A) which, when applied to a
%       harmonic field defined over the IP region, yields a vector of the
%       extrapolated field over the "external point" region (defined by the 3D
%       binary array maskEp).  
% 
%       Output argument M is a truncation (masking) operator for the EP, such that
%       M*x 'picks out' the EP from vector x.
%
%       Viz., to extrapolate a (3D but vectorized) reduced harmonic field (x) 
%       to the EP region defined by M: 
%
%       
%               = M' * A * x  
%
%
% =========================================================================
disp('########################') ;
disp('Extending harmonic field') ;
disp(['Expansion order: ' int2str(Parameters.expansionOrder)]) ;

gridSizeImg  = size( maskIp ) ;
nVoxelsImg = prod( gridSizeImg ) ;
nIp        = nnz( maskIp(:) ) ;
nEp        = nnz( maskEp(:) ) ;

%%
%%------- Map IP to EP
    
fout = fopen([Parameters.tmpSaveFldr 'maskEp.bin'], 'w') ;
fwrite(fout, uint8( maskEp(:) ) ) ;
fclose( fout ) ;

fout = fopen([Parameters.tmpSaveFldr 'maskIp.bin'], 'w') ;
fwrite(fout, uint8( maskIp(:) ) ) ;
fclose( fout ) ;

mapip2epArguments = [ ...
                  '  ' num2str( gridSizeImg ) ...
                  '  ' num2str( Parameters.radius ) ...
                  '  ' [Parameters.tmpSaveFldr 'maskIp.bin ' ] ...
                  '  ' [Parameters.tmpSaveFldr 'maskEp.bin ' ] ...
                  '  ' Parameters.tmpSaveFldr ] 

system([ 'mapip2ep ' mapip2epArguments]) ;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% LOAD FILES
% .../Tmp/displacementsX.bin
% .../Tmp/displacementsY.bin
% .../Tmp/displacementsZ.bin
% .../Tmp/rows.bin
% .../Tmp/columns.bin
% .../Tmp/indexIP.bin
Dx = fread( fopen( ...
	[ Parameters.tmpSaveFldr 'displacementsX.bin'], 'r'), inf, 'int32');

Dy = fread( fopen( ...
	[ Parameters.tmpSaveFldr 'displacementsY.bin'], 'r'), inf, 'int32' );

Dz = fread( fopen( ...
	[ Parameters.tmpSaveFldr 'displacementsZ.bin'], 'r'), inf, 'int32' );

rows = fread( fopen( ...
	[ Parameters.tmpSaveFldr 'rows.bin'], 'r'), inf, 'int32' );

columns = fread( fopen( ...
	[ Parameters.tmpSaveFldr 'columns.bin'], 'r'), inf, 'int32' );

indexEp = fread( fopen( ...
	[ Parameters.tmpSaveFldr 'indexEP.bin'], 'r'), inf, 'int32' ) ;

%% Adjust C zero-based indexing to MATLAB style:
%indexIP = indexIP + 1 ; 
rows    = rows + 1 ;
columns = columns + 1 ;
indexEp = indexEp + 1 ;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
disp('Deleting tmp files')
delete([ Parameters.tmpSaveFldr 'maskIp.bin']  ) ;
delete([ Parameters.tmpSaveFldr 'maskEp.bin']  ) ;
delete([ Parameters.tmpSaveFldr 'displacementsX.bin']);
delete([ Parameters.tmpSaveFldr 'displacementsY.bin']);
delete([ Parameters.tmpSaveFldr 'displacementsZ.bin']);
delete([ Parameters.tmpSaveFldr 'rows.bin']);
delete([ Parameters.tmpSaveFldr 'columns.bin']);
delete([ Parameters.tmpSaveFldr 'indexEP.bin']);

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% Create matrix operator A

nRows = nEp ;

Dx = sparse( rows, columns, Dx, nRows, nVoxelsImg ) ;
Dy = sparse( rows, columns, Dy, nRows, nVoxelsImg ) ;
Dz = sparse( rows, columns, Dz, nRows, nVoxelsImg ) ;

% (ideally: nEpRecoverable = nRows = nEp, but this depends on the chosen radius of expansion)

nEpRecoverable = length( unique( indexEp  ) ) ; 

disp(['Recovering ' num2str(nEpRecoverable) ' of total ' num2str(nEp) ...
    ' ( ' num2str(100*nEpRecoverable/nEp) '% ) ']);

if nEpRecoverable ~= nEp
    disp('Some EP remain out of reach with the current radius of expansion!') ;
end

clear rows columns

I  = ( abs(Dx) + abs(Dy) + abs(Dz) ) ~= 0 ;

% IP to EP distance:
W = ( ((Parameters.voxelSize(1)*Dx).^2 ) + ...
      ((Parameters.voxelSize(2)*Dy).^2 ) + ...
      ((Parameters.voxelSize(3)*Dz).^2 ) ) .^ 0.5 ;

% inverse as weights: W
W(I) = W(I) .^-1 ;

if Parameters.expansionOrder == 0
    disp('clearing 2') 
    clear Dx Dy Dz
end

% normalize:
W = spdiags( 1./sum(W, 2), 0, nRows, nRows ) * W ;

% matrix (Taylor expansion) operator:	
A = W .* I ;

if (Parameters.expansionOrder > 0)

	[Gx, Gy, Gz] = createdifferenceoperators( gridSizeImg, [1 1 1], 1 ) ;

	A  = A + [W .* Dx, W .* Dy, W .*Dz]*[Gx; Gy; Gz] ;

	if (Parameters.expansionOrder > 1)		
		Wx = (W .* Dx)/2 ;
		A  = A + [Wx .* Dx, Wx .* Dy, Wx .*Dz] * [Gx*Gx; Gy*Gx; Gz*Gx] ; 
		
		Wy = (W .* Dy)/2 ;
		A  = A + [Wy .* Dx, Wy .* Dy, Wy .*Dz] * [Gx*Gy; Gy*Gy; Gz*Gy] ; 
		
		Wz = (W .* Dz)/2 ;
		A  = A + [Wz .* Dx, Wz .* Dy, Wz .*Dz] * [Gx*Gz; Gy*Gz; Gz*Gz] ; 
	end
end

% truncation (masking) operator (e.g. M*x, 'picks out' the EP from vector x)
M = sparse( find(sum(I, 2) ~= 0 ), unique(indexEp), ones([nEpRecoverable 1]), nRows, nVoxelsImg ) ;

%M = sparse( 1:nRows, find( maskEp(:) ), ones([nRows 1]), nRows, nVoxelsImg ) ;





end

