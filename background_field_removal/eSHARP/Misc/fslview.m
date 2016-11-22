function [] = fslview( dataArray, Options, varargin )
%FSLVIEW view 3d array in fslview
%
%   Syntax
%
%   fslview( dataArray )
%   fslview( dataArray1, Options )
%   fslview( dataArray1, Options, dataArray2, dataArray3, ... )
%
%   Description
%
%   creates temporary nifti conversion of dataArray and passes it to fslview
%
%   
%       Options (object) to fslview input 
%       (see manual 'man fslview')
%       
%   The following Option-fields are supported
%
%   l (color scheme: Full list avaible in FSLView GUI)  
%       (e.g. one of: Greyscale, Red-Yellow, Blue-Lightblue,
%              Red, Green, Blue, Yellow, Pink, Hot, Cool, Copper, ...)
%       default: Options.l = 'Greyscale'
%
%   b (brightness)
%           if scalar S: brightness {-Abs|S|,+Abs|S|} 
%           if 2 component vector: brightness {S(1),S(2)} 
%           default: Options.b = [min(dataArray), max(dataArray)] 
%
%
% Ryan Topfer 2012
% topfer@ualberta.ca

% TODO
% variable brightness and colorschme for each overlay

img   = make_nii( double( dataArray ) ) ;

tmpImgNames = {'view3dtmp1.nii'} ;

save_nii(img, char(tmpImgNames(1))) ;

% call2Fslview = 'view3dtmp.nii ' ;

if nargin == 1 || isempty(Options)
    Options.dummy = [] ; %using defaults
end


%% brightness levels

if ~myisfield( Options, 'b' ) || isempty(Options.b)
    Options.b = [min(dataArray(:)) max(dataArray(:))] ;
    
elseif size( Options.b, 2) == 1 ;
    Options.b(2) = abs(Options.b) ;
    Options.b(1) = -abs(Options.b(1)) ;
end

viewOptions = [' -b "' num2str(Options.b(1)) ',' num2str(Options.b(2)) '"'];


%% color scheme

if myisfield( Options, 'l' )
    if isempty(Options.l)
        Options.l = 'Greyscale' ;
    end
    viewOptions = [ viewOptions ' -l ' Options.l ] ;
end

%%

call2Fslview = strcat( tmpImgNames(1), viewOptions ) ;

%% Additional input images?
numOverlays = nargin - 2 ;

for overlay = 1 : numOverlays
    
    img = make_nii( double( varargin{overlay} ) ) ;
    
    tmpImgNames{(overlay + 1)} = ['view3dtmp' int2str( overlay + 1 ) '.nii'] ;
    
    save_nii(img, char(tmpImgNames(overlay + 1)) ) ;

    call2Fslview = strcat( call2Fslview, {' '}, tmpImgNames{(overlay + 1)}, viewOptions) ;
end
    
call2Fslview = ['fslview ', char(call2Fslview), ' & '] ;

system( call2Fslview ) ;
pause(5)

call2rm = char( tmpImgNames ) ;
call2rm(:, end + 1) = ' ' ;

tmp = [] ;
for k = 1 : size( call2rm, 1 ) 
    tmp = [ tmp call2rm(k, :) ] ;
end

call2rm = tmp ;

system( ['rm ' call2rm] ) ;

end
