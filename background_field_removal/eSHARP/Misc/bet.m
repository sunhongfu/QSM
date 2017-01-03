function[ mask, f, skull ] = bet(mag, Options)
%BET - brain extraction tool
%
%   calls on BET of FSL package 
%       
%   Syntax
%
%   BET(Mag)
%   BET(Mag, Options)
%
%   Description
%       (see BET man page)
%
%   Mask            = BET(Mag)                   
%       returns binary brain mask based on (single-echo) 3D magnitude input 
%
%   Mask, f         = BET(Mag, Options)
%       also returns additional fractional intensity threshold 
%
%   Mask, f, skull  = BET(Mag, Options)
%       also returns binary mask of skull (might not be useful)
%
%   The following Option-fields are supported
%
%       voxelSize
%               default: [1 1 1] (isotropic)
%
%       f 
%           "fractional intensity threshold"
%               default: 0.5
% 
%   2013 topfer@ualberta.ca 


%% constants
DEFAULT_VOXELSIZE = [1 1 1] ;
DEFAULT_F         = 0.5 ;

%% check inputs

if nargin < 1 || isempty(mag)
    error('Function requires at least 1 input.')
end

if nargin < 2 || isempty(Options)
    disp('Default parameters will be used')
    disp('assuming isotropic resolution')
    Options.dummy = [] ;
end

if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end

if ~myisfield( Options, 'f' ) || isempty(Options.f)
    Options.f = DEFAULT_F ;
end


tmpFldr = '~/tmpBETFldr/' ;
system(['mkdir ' tmpFldr]) ;

save_nii( make_nii( mag, Options.voxelSize ), [tmpFldr 'mag.nii'] ) ;


betArgStr = [tmpFldr 'masked.nii -m -f ' num2str(Options.f) ] ;

if nargout == 3
    betArgStr = strcat( betArgStr, ' -s ' ) ;
end
    
system(['bet ' tmpFldr 'mag.nii ' betArgStr] ) ;


system(['gunzip ' tmpFldr 'masked_mask.nii.gz -d']) ;
mask = load_nii([tmpFldr 'masked_mask.nii']) ;
mask = single( mask.img ) ;

if nargout == 3
    system(['gunzip ' tmpFldr 'masked_skull.nii.gz -d']) ;
    skull = load_nii([tmpFldr 'masked_skull.nii']) ;
    skull = single( skull.img ) ;
end


system(['rm -r ' tmpFldr]) ;

f = Options.f ;


end