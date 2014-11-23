function[unwrappedPhase]=lapunwrap(rawPhase, Options)
%LAPUNWRAP phase unwrapping using Laplacian operator
%
%   LAPUNWRAP returns unwrapped phase array object  
%
%   Syntax
%
%   LAPUNWRAP(rawPhase)
%   LAPUNWRAP(rawPhase,Options)
%
%   Description
%
%   UnwrappedPhaseObj = LAPUNWRAP(rawPhase,[Options]) returns object with
%   identical fields to the input data object.
%
%
%   The following Option-fields are supported
%       voxelSize
%           default: [1 1 1] isotropic
%
% Based on Schofield & Zhu "Fast phase unwrapping algorithm for interferometric applications" Optics Letters 2003 

%% constants
DEFAULT_VOXELSIZE = [1 1 1];


%% check inputs

if nargin < 2 || isempty(Options)
    disp('assuming isotropic resolution')
    Options.dummy = [] ;
end

if  ~myisfield( Options, 'voxelSize' ) || isempty(Options.voxelSize)
    Options.voxelSize = DEFAULT_VOXELSIZE ;
end

% forloop to unwrap input phase in more than 3D formats
[np, nv, nv2, numEcho, numRx] = size(rawPhase);
unwrappedPhase = zeros(size(rawPhase));

GDV   = [np,nv,nv2];
tmp   = zeros( 2*GDV );

k2    = sum( makekspace( 2*GDV, Options ).^2, 4) ;
k2    = fftshift( k2 ) ;

for i = 1:numEcho*numRx
    %% initialize
    phase = squeeze(rawPhase(:,:,:,i)) ;
    
    %% mirroring
    tmp(1:GDV(1), 1:GDV(2), 1:GDV(3))                            = phase ;
    tmp(1:GDV(1), GDV(2)+1:2*GDV(2), 1:GDV(3))                   = flipdim(phase,2);
    tmp(1:GDV(1), 1:GDV(2), GDV(3)+1:2*GDV(3))                   = flipdim(phase,3);
    tmp(1:GDV(1), GDV(2)+1:2*GDV(2), GDV(3)+1:2*GDV(3))          = flipdim(flipdim(phase,2),3);
    tmp(GDV(1)+1:2*GDV(1), GDV(2)+1:2*GDV(2), GDV(3)+1:2*GDV(3)) = flipdim(flipdim(flipdim(phase,1),2),3);
    tmp(GDV(1)+1:2*GDV(1), GDV(2)+1:2*GDV(2), 1:GDV(3))          = flipdim(flipdim(phase,1),2);
    tmp(GDV(1)+1:2*GDV(1), 1:GDV(2), GDV(3)+1:2*GDV(3))          = flipdim(flipdim(phase,1),3);
    tmp(GDV(1)+1:2*GDV(1), 1:GDV(2), 1:GDV(3))                   = flipdim(phase,1);
    
    %% unwrap
    phase = tmp; 
    clear tmp

    phase = fftn(cos(phase).*ifftn(k2.*fftn(sin(phase)))-sin(phase).*ifftn(k2.*fftn(cos(phase))))./(k2+eps);
    
    phase( abs(phase) == Inf ) = 0;
    
    phase = ifftn(phase);
    
    %% output
    unwrappedPhase(:,:,:,i) = real(phase(1:GDV(1),1:GDV(2),1:GDV(3)));
end

