classdef cls_tv_new 
%CLS_TV Class for total variation
    
    properties
        adjoint = 0
        mask
    end
    
    methods
        function this = cls_tv_new(mask)
            this.mask   = mask;
        end

        function obj = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
        end
        
        function product = mtimes(obj,b)
            if obj.adjoint
                product = obj.mask.*invD((1-repmat(obj.mask,[1 1 1 3])).*b);
            else
                product = D(obj.mask.*b);
            end
        end
    end
    

        
end


function res = invD(y)
    res = adjDx(y(:,:,:,1)) + adjDy(y(:,:,:,2)) + adjDz(y(:,:,:,3));
end

function res = D(x)
    Dx = x([2:end,end],:,:) - x;
    Dy = x(:,[2:end,end],:) - x;
    Dz = x(:,:,[2:end,end]) - x;
    res = cat(4,Dx,Dy,Dz);
end

function res = adjDx(x)
    res = x([1,1:end-1],:,:) - x;
    res(1,:,:) = -x(1,:,:);
    res(end,:,:) = x(end-1,:,:);
end

function res = adjDy(x)
    res = x(:,[1,1:end-1],:) - x;
    res(:,1,:) = -x(:,1,:);
    res(:,end,:) = x(:,end-1,:);
end

function res = adjDz(x)
    res = x(:,:,[1,1:end-1]) - x;
    res(:,:,1) = -x(:,:,1);
    res(:,:,end) = x(:,:,end-1);
end

