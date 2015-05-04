function[dataArray] = makeodd(dataArray, option, finalGridDimensionVector)
%MAKEODD returns dataArray with odd dimensions
%   
%   Syntax
%
%   A = MAKEODD(A)
%   A = MAKEODD(A, option)
%   A = MAKEODD(A, option, gridDimensionVector) ;
%
%   returns array A with odd dimensions.
%   
%   option
%
%       if 'isPadding', even dimensions are padded with a slice of zeros 
%           (default: even dimensions have a slice shaved off)
%
%       if 'isUndoing', odd dimensions are padded with a slice of zeros
%           (default: make ODD!)
%
%       gridDimensionVector
%           specifies final array size
%           (default: minimum change from original size such that 'option'
%           is satisfied)
%

if nargin < 2
    option = 'isShaving' ;
    isPadding = false ;
end

if strcmp(option, 'isUndoing')
    isPadding = true ;
end


%%

gridDimensionVector = size( dataArray ) ;

isOdd               = mod( gridDimensionVector, 2 ) ;
isEven              = double(~isOdd) ;

%%

if nargin == 3
    numSlices           = abs( gridDimensionVector - finalGridDimensionVector ) ;
else
    
    switch option
        case 'isShaving'
            numSlices = isEven ;
        case 'isUndoing'
            numSlices = isOdd ;
    end
end
    
    
if isPadding
    
    dataArray = padarray(dataArray, numSlices, 'post') ;
     
else
    
    dataArray = dataArray( 1:end-numSlices(1), 1:end-numSlices(2),1:end-numSlices(3)) ;
        
end

end




