function IM3D=unDCshift(IM3D,Pars)

% Removes DC shift of a raw image matrix by subtracting the mean

% Produced by Sandra Meyers/Amir Eissa
% June 19/09

% Input: complex 3D matrix and parameters variable

Mean=squeeze(mean(IM3D(1:round(Pars.np/20),:,:)));
DCshiftmat=repmat(Mean,[1 1 Pars.np/2]);
clear Mean
DCshiftmat=shiftdim(DCshiftmat,2);
IM3D=IM3D-DCshiftmat; 