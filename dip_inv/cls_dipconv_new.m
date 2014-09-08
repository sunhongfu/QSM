% classdef cls_dipconv_new
% %CLS_DIPCONV Class for unit dipole kernel convolution
    
%     properties
%         transp = 0
%         imsize
%         ker  % unit dipole kernel in k-space 'D'
%         mask
%         res
%     end
    
%     methods
%         function this = cls_dipconv_new(imsize,ker,mask,res)
%             this.imsize = imsize;
%             this.ker = ker;
%             this.mask = mask;
%             this.res = res;
%         end
        
%         function obj = ctranspose(obj)
%             obj.transp = xor(obj.transp,1);
%         end
        
%         function y = mtimes(obj,x)
%             x = reshape(x,obj.imsize);
%             if obj.transp
%                 y = (1-obj.mask).*ifftn(obj.ker.*fftn(x));
%             else
%                 y = ifftn(obj.ker.*fftn(obj.mask.*obj.res+(1-obj.mask).*x));
%             end

            
%         end
%     end
    
% end

classdef cls_dipconv_new
%CLS_DIPCONV_new Class for unit dipole kernel convolution
    
    properties
        transp = 0
        imsize
        ker  % unit dipole kernel in k-space 'D'
        mask_hemo
    end
    
    methods
        function this = cls_dipconv_new(imsize,ker,mask_hemo)
            this.imsize = imsize;
            this.ker = ker;
            this.mask_hemo = mask_hemo;
        end
        
        function obj = ctranspose(obj)
            obj.transp = xor(obj.transp,1);
        end
        
        function y = mtimes(obj,x)
            x = reshape(x,obj.imsize);
            if obj.transp
                y = obj.mask_hemo.*ifftn(obj.ker.*fftn(x));
            else
                y = ifftn(obj.ker.*fftn(obj.mask_hemo.*x));
            end
        end
    end
    
end

