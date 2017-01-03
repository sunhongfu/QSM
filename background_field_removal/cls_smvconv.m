classdef cls_smvconv
%CLS_SMVCONV Class for SMV convolution
    
    properties
        transp = 0 % transpose flag
        mask       % region of interest / binary mask
        csh        % circular shift matrix
        ker        % deconvolution kernel in k-space 'C'
        imsize
    end
    
    methods
        function this = cls_smvconv(imsize,ker,csh,mask)
            this.imsize = imsize;
            this.ker    = ker;
            this.csh    = csh;
            this.mask   = mask;            
        end
        
        function obj = ctranspose(obj)
            obj.transp = xor(obj.transp,1);
        end
        
        function y = mtimes(obj,x) % x,y are column vectors
            x = reshape(x,obj.imsize);
            if obj.transp
                y = ifftn(conj(obj.ker).*fftn(circshift(obj.mask.*x,obj.csh)));
                y = y(:);
            else
                y = obj.mask.*circshift(ifftn(obj.ker.*fftn(x)),-obj.csh);
                y = y(:);
            end
        end
    end
    
end

