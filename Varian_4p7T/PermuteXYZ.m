function Output = PermuteXYZ(Input,Pars)

% function Output = PermuteXYZ(Input,Pars)
% This function permutes a 3D data matrix (even if acquired as 2D set) into
% The order of (X,Y,Z) assuming that no two acquisition size numbers are 
% identical

ReadOutNumber=Pars.np/2; % not the right case if ACQ is done efficiently while 
                    % partial read out is ON ... it is fine for now
PhaseYNumber=Pars.nv;
PhaseZNumber=max(Pars.nv2,Pars.NS);

S=size(Input);

OrderX=find(S==ReadOutNumber);
OrderY=find(S==PhaseYNumber);
OrderZ=find(S==PhaseZNumber);

Output=permute(Input,[OrderX OrderY OrderZ ]);

return