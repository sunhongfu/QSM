function[imgOut]=resizekspace(imgIn, resizeFactor)
%RESIZEKSPACE
%
%   Syntax
%
%   Img = RESIZEKSPACE(Img, resizeFactor)
%
%   Description
%
%   Img = RESIZEKSPACE(Img, resizeFactor)
%       
%       if resizeFactor < 1, RESIZEKSPACE extracts centre of Img k-space to return low
%   resolution image.
%       if resizeFactor > 1, 
%
%   Input
%       Img
%           complex image. First 3 dimensions should be spatial x,y, z.
%           Optional 4th and 5th dimensions can refer to echo and channel.
%
%       resizeFactor
%           1 or 2 or 3 component vector
%               e.g. = .5, returns image with half the resolution in each
%               dimension
%               e.g. = [1 2 1] returns interpolated image (2x original
%               resolution in 2nd dimension)
%
%   Note: expects each 3D image to possess EVEN dimensions
%
%
% TODO: correct output scaling

%%check inputs
if nargin < 2 || isempty(imgIn) || isempty(resizeFactor)
    error('Function requires at least 2 input arguments.')
end

if numel(resizeFactor) == 1
    disp('single element resize factor given')
    disp(['resizing all dimensions by factor ' num2str(resizeFactor)])   
    resizeFactor(2:3) = resizeFactor(1) ;
end

if numel(resizeFactor) == 2
    disp('2 element resize factor given')
    disp('resolution of 3rd dimension will not be changed')
    resizeFactor(3) = 1 ;
end

if any( resizeFactor < 1 )
    isCroppingKSpace  = true ;
else
    isCroppingKSpace  = false ;
end



gridDimensionVector = size( abs( imgIn ) ) ;

if numel(gridDimensionVector) > 3 ;
    numEcho = gridDimensionVector(4) ;
else
    numEcho = 1 ;
end

if numel(gridDimensionVector) > 4
    numRx   = gridDimensionVector(5) ;
else
    numRx = 1 ;
end

gridDimensionVector = gridDimensionVector(1:3) ;


newGridDimensionVector = round(gridDimensionVector .* resizeFactor) ;

imgOut = zeros( [newGridDimensionVector numEcho numRx] ) ;

if max( gridDimensionVector - newGridDimensionVector ) > 0

    disp('cropping k-space...')
    
    % even-dimension array...
    midPoint   = gridDimensionVector/2 ;
    
    offset      = midPoint - round(midPoint .* resizeFactor) + 1 ;
    offset(2,:) = midPoint + round(midPoint .* resizeFactor) ;
    
    for echo = 1 :numEcho
        for rx = 1 : numRx
            
            fImg = fftshift( fftc( imgIn(:,:,:, echo, rx) ) ) ;
            
            fImg = fImg(offset(1,1):offset(2,1), offset(1,2) : offset(2,2),offset(1,3) : offset(2,3));
            
            img  = ifftc( ifftshift( fImg ) ) ;

            imgOut(:,:,:,echo,rx) = img ;
            
        end
    end
    
end

if min( gridDimensionVector - newGridDimensionVector ) < 0
    disp('zero-padding in k-space...')
    
    if isCroppingKSpace
        imgIn = imgOut ;
    end
    
    padSize = round( (newGridDimensionVector - gridDimensionVector )/2 ) ;
    
    
    
    for echo = 1 :numEcho
        for rx = 1 : numRx
            
            fImg = fftshift( fftc( imgIn(:,:,:, echo, rx) ) ) ;
            
            fImg = padarray( fImg, padSize ) ;
            

            img  = ifftc( ifftshift( fImg ) ) ;
            

            imgOut(:,:,:,echo,rx) = prod(resizeFactor)*img ;
            
        end
    end
    
end
