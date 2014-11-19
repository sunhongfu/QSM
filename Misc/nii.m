function[] = nii(dataArray, name, Options)
%NII 
%
%   NII
%
%
%   Syntax
%
%   NII(A)
%   NII(A,name)
%   NII(A,name,Options)
%   .......................
%
%   Description
%   
%   
%
%   .......................
%   
%   Input 
%   
%   A 
%       data array being saved as nifti
%
%   name
%       file name base (e.g. 'BigSmiley' results in file 'BigSmiley.nii')
%
%   Options
%       object with possible fields:
%           .voxelSize
%               (default: [1 1 1])
%      
%           .isReorienting
%               (true [default]||false)
%               if true, array is flipped to have longest dimension (e.g. read) 
%               along the vertical
%   
%           .mediaSaveFldr
%               (default: .)
%               


DEFAULT_MEDIASAVEFLDR = './' ;
DEFAULT_ISREORIENTING = true ;
DEFAULT_NAME          = 'tmp' ;
DEFAULT_VOXELSIZE     = [1 1 1] ;


if nargin < 1 || isempty(dataArray) 
    error('Function requires at least 1 input argument (2d or 3d array).')
end

if nargin < 2 || isempty(name)
    name = DEFAULT_NAME ;
end

if nargin < 3 || isempty(Options)
    Options.dummy = [] ;
end

if  ~myisfield( Options, 'mediaSaveFldr' ) || isempty(Options.mediaSaveFldr)
    Options.mediaSaveFldr = DEFAULT_MEDIASAVEFLDR ;
end

if  ~myisfield( Options, 'isReorienting' ) || isempty(Options.isReorienting)
    Options.isReorienting = DEFAULT_ISREORIENTING ;
end

if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end




gridDimensionVector  = size( dataArray ) ;

if numel( gridDimensionVector ) < 4
    dataArray(:,:,:,1) = dataArray ;
end

if Options.isReorienting
    [tmp, readDimension] = max( gridDimensionVector ) ;
    
    if readDimension == 1
        dataArray = permute( dataArray, [2 readDimension 3 4] ) ;
        dataArray = flipdim(dataArray, 1) ;
    end
end





save_nii( make_nii(dataArray, Options.voxelSize), [ Options.mediaSaveFldr name '.nii']) ;

end