function [ unwrappedPhase ] = unwrap3d( rawPhase, mask, Options )
%UNWRAP3D 
%
%   UNWRAP3D calls a 3D path-based unwrapper (Abdul-Rahman et al 2007)
%
%	unwrappedPhase = UNWRAP3D( rawPhase, mask )
%	unwrappedPhase = UNWRAP3D( rawPhase, mask, Options )
%
%    .......................
%    
%    rawPhase: 3d (double) array of wrapped phase
%    mask:     3d (binary) array demarking region (of reliable raw phase
%               measurements) to be unwrapped
%
%    .......................
%
%   
%   The following Option-fields are supported
%
%        .path2UnwrappedPhase 
%               default: './unwrappedPhase.bin' 
%
%        .isSavingBinaries
%                   The following files are created in the dir of 
%                   path2UnwrappedPhase: 
%                      rawPhase.bin [float32], 
%                      unwrappedPhase.bin [float32],
%                      mask.bin [uint8]
%               default: false
%                   .bin files deleted. Assumption being the user simply wants
%                   the unwrapped phase returned as a matlab array.
%
%    .......................
% 	topfer@ualberta.ca	2015
%
path2HusseinExecutable = '~/Projects/General/bin/Hussein_3D_unwrapper_with_mask_backup' ;

DEFAULT_PATH2UNWRAPPEDPHASE   = [ './unwrappedPhase.bin' ] ;
DEFAULT_ISSAVINGBINARIES      = false ;


% check inputs

if nargin < 2 || isempty(rawPhase) || isempty(mask) 
    error('Function requires at least 2 input arguments.')
end

if nargin == 2
    Options.dummy = [] ;
end

if  ~myisfield( Options, 'path2UnwrappedPhase' ) || isempty(Options.path2UnwrappedPhase)
    Options.path2UnwrappedPhase = DEFAULT_PATH2UNWRAPPEDPHASE ;
end

if  ~myisfield( Options, 'isSavingBinaries' ) || isempty(Options.isSavingBinaries)
    Options.isSavingBinaries = DEFAULT_ISSAVINGBINARIES;
end



[dataSaveDirectory,~,~] = fileparts( Options.path2UnwrappedPhase ) ;
dataSaveDirectory = [dataSaveDirectory '/'] ;


% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

% permute so .bin images retain the same orientation when loaded 
% into imageJ as the original arrays
rawPhase = single( permute( rawPhase, [ 2 1 3] ) ) ;

% for the unwrapper, a voxel of value 255 is unwrapped, others ignored.
mask     = permute( 255*uint8(mask), [2 1 3] ) ; 

gridSize = size( rawPhase ) ;
nVoxels  = prod( gridSize ) ;

% =========================================================================
% write binaries 

fOut = fopen([dataSaveDirectory 'rawPhase.bin'], 'w') ;
fwrite( fOut, rawPhase(:), 'single' );
fclose( fOut ) ;

fOut = fopen([dataSaveDirectory 'mask.bin'], 'w') ;
fwrite( fOut, mask(:), 'uint8' );
fclose( fOut ) ;

%-------
% send command

unwrapCmd = [path2HusseinExecutable ...
             ' ' dataSaveDirectory 'rawPhase.bin' ...
             ' ' dataSaveDirectory 'mask.bin' ...
             ' ' dataSaveDirectory 'unwrappedPhase.bin' ...
             ' ' num2str( gridSize ) ] ;

system( unwrapCmd ) ; 

% =========================================================================
% read binary & cleanup

fIn = fopen([Options.path2UnwrappedPhase],'r') ;
unwrappedPhase = fread( fIn, nVoxels, 'single=>single' ) ;
fclose( fIn ) ;

unwrappedPhase = reshape( unwrappedPhase, gridSize )  ;
unwrappedPhase = (mask~=0).*unwrappedPhase ;
unwrappedPhase = permute( unwrappedPhase, [2 1 3] ) ;

if ~Options.isSavingBinaries
    delete( [dataSaveDirectory 'rawPhase.bin'], ...
            [dataSaveDirectory 'mask.bin'], ...
            [dataSaveDirectory 'unwrappedPhase.bin'] ) ;
end


end
