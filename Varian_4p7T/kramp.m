function IM3D=kramp(Pars,IM3D)

% If data is 3D, add a ramp to shift k-space to proper location
% Also permute matrix so dimensions are in the same order as 2D (for
% purpose of filtering/viewing) i.e. (x,z,y)

if Pars.nv2>1    
    slice=Pars.nv2; 
    ind = -Pars.nv2/2:Pars.nv2/2-1;
    pramp = exp(-1i.*(2*pi*(Pars.pss+0.5*Pars.lpe2/Pars.nv2)/Pars.lpe2).*ind);

    kramp=ones(Pars.nv2,Pars.np/2);
    for iq=1:Pars.nv2
    kramp(iq,:)=kramp(iq,:).*pramp(iq);
    end
    Pars.nv2=slice;
    for i=1:Pars.nv
       temp(:,:)=squeeze(IM3D(:,i,:));
       if size(temp)==size(kramp)
           temp=temp.*kramp;
       else
           temp=temp.*kramp';
       end

       IM3D(:,i,:)=temp;
       
       temp=temp*0;

    end
    clear temp*

end
