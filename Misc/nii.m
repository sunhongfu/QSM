function img = nii(img, Options)
%NII 
%
%   NII
%
%
%   Syntax
%
%   NII(A)
%   NII(A,Options)
%   .......................
%
%   Description
%   
%   NII calls the make_nii() and save_nii() functions of the NifTI toolbox.
%   However the toolbox seems to reverse read/phase-encode (column/row) order
%   compared with that of the original image (DICOM). NII then first swaps 
%   the rows & columns such that there is no apparent rotation after writing to
%   file with save_nii().
%
%   .......................
%   
%   Input 
%   
%   A 
%       data array being saved as nifti
%
%   Options
%       object with possible fields:
%           .voxelSize
%               (default: [1 1 1])
%      
%           .filename
%               (without or without .nii extension)
%               (default: './tmp.nii')
%  
%  TODO DOCUMENTATION         

DEFAULT_VOXELSIZE = [1 1 1] ;
DEFAULT_FILENAME  = './tmp.nii' ;

if nargin < 1 || isempty(img) 
    error('Function requires at least 1 input argument (2d or 3d array).')
end

if nargin < 2 || isempty(Options)
    Options.dummy = [] ;
end

isConverting2Nii = true ;


tmp = whos('img') ;

if strcmp( tmp.class, 'char' )
    isConverting2Nii = false ; % a nifti image is to be loaded and returned as a matlab array
end

if isConverting2Nii
    if  ~myisfield( Options, 'filename' ) || isempty(Options.filename)
        Options.filename = DEFAULT_FILENAME ;
    else
        if length( Options.filename < 4) || ( ~strcmp( Options.filename(end-3:end), '.nii') )
            Options.filename = [Options.filename '.nii' ] ;
        end
    end

    if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
        Options.voxelSize = DEFAULT_VOXELSIZE ;
    end

    img = double(img) ;

    % flip rows & columns for make_nii
    img = permute( img, [2 1 3 4] ) ;

    save_nii( make_nii(img, Options.voxelSize), Options.filename) ;

else

    img = load_nii( img ) ;
    img = double( img.img ) ;

    % flip rows & columns for load_nii
    img = permute( img, [2 1 3 4] ) ;


end
