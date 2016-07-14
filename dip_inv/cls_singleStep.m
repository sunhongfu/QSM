classdef cls_singleStep
%CLS_DIPCONV Class for unit dipole kernel convolution
    
    properties
        transp = 0 % transpose flag
        mask
        DKER  
        D    % unit dipole kernel in k-space 'D'
        csh
    end
    
    methods
        function this = cls_singleStep(mask,DKER,D,csh)
            this.mask = mask;
            this.DKER = DKER;
            this.D    = D;
            this.csh  = csh;
        end
        
        function obj = ctranspose(obj)
            obj.transp = xor(obj.transp,1);
        end
        
        function y = mtimes(obj,x) % x is column vector
            x = reshape(x,size(obj.mask));
            if obj.transp
                y = ifftn(conj(obj.D).*conj(obj.DKER).*fftn(circshift(obj.mask.*x,obj.csh)));
                % y = y(:);
            else
                y = obj.mask.*circshift(ifftn(obj.DKER.*obj.D.*fftn(x)),-obj.csh);
                % y = y(:);
            end
        end
    end
    
end

