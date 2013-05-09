function[dataArray] = makeodd(dataArray, option)

%MAKEODD returns dataArray with odd dimensions
%   
%   Syntax
%
%   A = MAKEODD(A)
%   A = MAKEODD(A, option)
%
%   returns array A with odd dimensions
%   
%   option
%
%       if 'isPadding', even dimensions are padded with a slice of zeros 
%           (default: even dimensions have a slice shaved off)
%
%       if 'isUndoing', odd dimensions are padded with a slice of zeros
%           (default: make ODD!)


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


if isPadding
        
        if strcmp(option,'isUndoing')
            %pad array to even dimensions
            dataArray = padarray(dataArray, isOdd, 'post') ;
        else
            %pad array to odd dimensions
            dataArray = padarray(dataArray, isEven, 'post') ;
        end
    
        
else
        if isEven(1)
            dataArray = dataArray( 1:end-1, :,:) ;
        end
        if isEven(2)
            dataArray = dataArray( :,1:end-1, :) ;
        end
        if isEven(3)
            dataArray = dataArray( :,:,1:end-1) ;
        end
        
        
end

end




