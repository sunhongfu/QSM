classdef cls_dipconv_polyfit
    
    properties
        transp = 0 % transpose flag
        imsize
        ker
        sus_mask
        polyn % the polynomials
    end
    
    methods
        function this = cls_dipconv_polyfit(imsize,ker,sus_mask,polyn)
            this.imsize   = imsize;
            this.ker      = ker;
            this.sus_mask = sus_mask;   
            this.polyn    = polyn;        
        end
        
        function obj = ctranspose(obj)
            obj.transp = xor(obj.transp,1);
        end
        
        function y = mtimes(obj,x) % x,y are column vectors
            if obj.transp
                top = obj.sus_mask.*ifftn(obj.ker.*fftn(x));
                bottom = (obj.polyn)'*x(:);
                y = [top(:); bottom];
            else
                chi = reshape(x(1:prod(obj.imsize)),obj.imsize);
                coeff = x(prod(obj.imsize)+1:end);
                y = ifftn(obj.ker.*fftn(obj.sus_mask.*chi)) + reshape(obj.polyn*coeff,obj.imsize);
            end
        end
    end
    
end

