classdef cls_dipconv
%CLS_DIPCONV Class for unit dipole kernel convolution
    
    properties
        imsize
        ker  % unit dipole kernel in k-space 'D'
    end
    
    methods
        function this = cls_dipconv(imsize,ker)
            this.imsize = imsize;
            this.ker = ker;
        end
        
        function obj = ctranspose(obj)
            % empty function
            % (fDF)' = fDF
        end
        
        function y = mtimes(obj,x)
            x = reshape(x,obj.imsize);
            y = ifftn(obj.ker.*fftn(x));
        end
    end
    
end

