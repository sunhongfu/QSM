function [ unwrappedPhase ] = prelude( rawPhase, mag, Options )
%PRELUDE 
%
%   PRELUDE calls FSL Prelude
%
%	unwrappedPhase = PRELUDE( rawPhase, mag )
%	unwrappedPhase = PRELUDE( rawPhase, mag, Options )
%
%    .......................
%
%    rawPhase: 3d array of the wrapped phase values
%    mag:      3d array of the corresponding magnitude data
%
%    .......................
%   
%   The following Option-fields are supported
%
%        .voxelSize
%               default: [1 1 1]
%
%        .path2UnwrappedPhase 
%               default: 'unwrappedPhase.nii' 
%                       temporarily created in the same directory as the raw
%                       phase input, however it is subsequently deleted
%                       (assumption being the user wants the unwrapped phase
%                       returned merely as a matlab array).
%
%        .mask
%
%
%        .isUnwrappingIn2D
%                        default: true
%
%        .isSavingNiftis
%                        default: false
%
%    .......................
% 	topfer@ualberta.ca	2014

	DEFAULT_VOXELSIZE             = [1 1 1] ;
	DEFAULT_PATH2UNWRAPPEDPHASE   = [ './unwrappedPhase.nii' ] ;
	DEFAULT_ISUNWRAPPINGIN2D      = true ;
	DEFAULT_ISSAVINGNIFTIS        = false ;

	% check inputs

	if nargin < 2 || isempty(rawPhase) || isempty(mag) 
	    error('Function requires at least 2 input arguments.')
	end

	if nargin < 3 || isempty(Options)
	    disp('Default parameters will be used')
	    Options.dummy = [] ;
	end

	if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
	    Options.voxelSize = DEFAULT_VOXELSIZE ;
	end

	if  ~myisfield( Options, 'path2UnwrappedPhase' ) || isempty(Options.path2UnwrappedPhase)
		isSavingUnwrappedPhase = false ;
	    Options.path2UnwrappedPhase = DEFAULT_PATH2UNWRAPPEDPHASE ;
	end

	if  ~myisfield( Options, 'isUnwrappingIn2D' ) || isempty(Options.isUnwrappingIn2D)
	    Options.isUnwrappingIn2D = DEFAULT_ISUNWRAPPINGIN2D;
	end

	if  ~myisfield( Options, 'isSavingNiftis' ) || isempty(Options.isSavingNiftis)
	    Options.isSavingNiftis = DEFAULT_ISSAVINGNIFTIS;
	end

	[dataSaveDirectory,~,~] = fileparts( Options.path2UnwrappedPhase ) ;
	dataSaveDirectory = [dataSaveDirectory '/'] ;

	NiiOptions.voxelSize = Options.voxelSize ;
	
	NiiOptions.filename  = [dataSaveDirectory 'rawPhase'] ;
	nii( rawPhase, NiiOptions ) ;
	
	NiiOptions.filename  = [dataSaveDirectory 'mag'] ;
	nii( mag, NiiOptions ) ;

	Options.etc = ' ' ;

	if Options.isUnwrappingIn2D
		Options.etc   = ' -s ' ;
	end

	if myisfield( Options, 'mask' )
		NiiOptions.filename  = [dataSaveDirectory 'mask'] ;
		nii( Options.mask, NiiOptions ) ;

		Options.etc   = [ Options.etc ' -m ' dataSaveDirectory 'mask' ] ;
	end

	unwrapCommand = ['prelude -p ' dataSaveDirectory 'rawPhase' ...
					 		' -a ' dataSaveDirectory 'mag' ...
					 		' -o ' Options.path2UnwrappedPhase ...
					 		' ' Options.etc] ; 

	system(unwrapCommand) ;

	system(['gunzip ' Options.path2UnwrappedPhase '.gz -df']) ;

	unwrappedPhase = nii( Options.path2UnwrappedPhase ) ;

	if ~Options.isSavingNiftis
		delete( [dataSaveDirectory 'rawPhase.nii'], [dataSaveDirectory 'mag.nii'], ...
				[dataSaveDirectory 'unwrappedPhase.nii'] ) ;
                
                if myisfield( Options, 'mask' )	
                        delete( [dataSaveDirectory 'mask.nii'] ) ;
                end
    end
end
