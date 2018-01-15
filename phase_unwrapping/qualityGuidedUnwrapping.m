%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Véronique Fortier 2017 - MIT License
%
% Quality guided unwrapping                                               
%
% Inputs:                                                                 
% -Phase image (im_phase)                                                 
% -Object mask (im_mask) - to reduce unwrapping time                                               
% -Quality threshold (qualityCutoff): This threshold has to be in the form
% of 100/(percent threshold). Default is 3.5 (equivalent to 29%)          
%
% Outputs:                                                                
% -Unwrapped result (im_unwrapped)                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ im_unwrapped ] = qualityGuidedUnwrapping(im_phase, im_mask,qualityCutoff )

im_unwrapped=zeros(size(im_phase));               % starting matrix for unwrapped phase
adjoin=zeros(size(im_phase));                     % starting matrix for adjoin matrix
unwrapped_binary=zeros(size(im_phase));           % Binary image to mark unwrapped pixels
matrixSize=size(im_mask);

if nargin<3
    qualityCutoff=3.5;
end


%% Calculate second difference quality map

im_phase_hor_right=zeros(matrixSize);
im_phase_hor_right(2:(end-1),2:(end-1),2:(end-1))=im_phase((2:(end-1))-1,2:(end-1),2:(end-1));
im_phase_hor_left=zeros(matrixSize);
im_phase_hor_left(2:(end-1),2:(end-1),2:(end-1))=im_phase((2:(end-1))+1,2:(end-1),2:(end-1));

im_phase_ver_top=zeros(matrixSize);
im_phase_ver_top(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),(2:(end-1))+1,2:(end-1));
im_phase_ver_bot=zeros(matrixSize);
im_phase_ver_bot(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),(2:(end-1))-1,2:(end-1));

im_phase_norm_back=zeros(matrixSize);
im_phase_norm_back(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),2:(end-1),(2:(end-1))-1);
im_phase_norm_front=zeros(matrixSize);
im_phase_norm_front(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),2:(end-1),(2:(end-1))+1);

H=wrapToPi(im_phase_hor_right-im_phase)-wrapToPi(im_phase-im_phase_hor_left);
V=wrapToPi(im_phase_ver_bot-im_phase)-wrapToPi(im_phase-im_phase_ver_top);
N=wrapToPi(im_phase_norm_back-im_phase)-wrapToPi(im_phase-im_phase_norm_front);

SD=(sqrt(H.^2+V.^2+N.^2));


%% Identify the starting seed point on the phase quality map (central region of the central slice)
im_phase_quality=SD.*im_mask;
SE=strel('square',5);
im_phase_qualityMask=im_phase_quality(:,:,round(matrixSize(3)/2)).*imerode(im_mask(:,:,round(matrixSize(3)/2)),SE);
im_phase_qualityMask=uint16(im_phase_qualityMask*100);
im_phase_qualityMask(find(im_phase_qualityMask==0))=1000;   % Set region outside of im_mask to an arbitrary very low quality 

[a, loc]=min(im_phase_qualityMask(:));  % The minimum is the highest quality point
[xpoint,ypoint]=ind2sub(size(im_phase_qualityMask),loc);
zpoint=round(matrixSize(3)/2);
 
 
%% Define the seed point as unwrapped
colref=round(ypoint);
rowref=round(xpoint);
im_unwrapped(rowref,colref,zpoint)=im_phase(rowref,colref,zpoint);                        
unwrapped_binary(rowref,colref,zpoint)=1;


%% Add the adjoining voxels of the seed point to the adjoin matrix
if im_mask(rowref-1, colref, zpoint)==1 adjoin(rowref-1, colref, zpoint)=1; end       
if im_mask(rowref+1, colref, zpoint)==1 adjoin(rowref+1, colref, zpoint)=1; end
if im_mask(rowref, colref-1, zpoint)==1 adjoin(rowref, colref-1, zpoint)=1; end
if im_mask(rowref, colref+1, zpoint)==1 adjoin(rowref, colref+1, zpoint)=1; end
if im_mask(rowref, colref, zpoint-1)==1 adjoin(rowref, colref, zpoint-1)=1; end
if im_mask(rowref, colref, zpoint+1)==1 adjoin(rowref, colref, zpoint+1)=1; end


%% Quality-guided unwrapping
im_unwrapped=GuidedFloodFill_3D(im_phase, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask,qualityCutoff);   

end

