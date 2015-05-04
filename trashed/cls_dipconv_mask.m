classdef cls_dipconv_mask
%CLS_DIPCONV Class for unit dipole kernel convolution
    
    properties
        imsize
        ker  % unit dipole kernel in k-space 'D'
        mask
    end
    
    methods
        function this = cls_dipconv_mask(imsize,ker,mask)
            this.imsize = imsize;
            this.ker = ker;
            this.mask = mask;
        end
        
        function obj = ctranspose(obj)
            % empty function
            % (fDF)' = fDF
        end
        
        function y = mtimes(obj,x)
            x = reshape(x,obj.imsize);
            y = obj.mask.*ifftn(obj.ker.*fftn(obj.mask.*x));
        end
    end
    
end

