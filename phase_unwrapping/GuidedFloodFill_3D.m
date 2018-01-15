%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Véronique Fortier 2017 - MIT License
%
% Quality guided unwrapping using a two sections and stack chain guiding
% strategy
%
% Input: 
% -Wrapped phase image (IM_phase) 
% -Unwrapped image (IM_unwrapped) - initially includes the seed point only
% -List of unwrapped pixels (unwrapped_binary) - initially includes the seed point only
% -Quality map (qualityMap)
% -An adjoining matrix (adjoin) - initially includes the pixels around the seed point only
% -An object mask (IM_mask) - To reduce unwrapping time
% -Quality threshold (qualityCutoff) - Default is 3.5
%
% Output:
% -Unwrapped result (IM_unwrapped)  
%
% Based on: GuidedFloodFill.m (2D unwrapping) - Created by B.S. Spottiswoode on 11/11/2004 
% (https://www.mathworks.com/matlabcentral/fileexchange/29499-qualityguidedunwrap2d-r1)
%
% Extended to 3D and accelerated with quality threshold by Veronique 
% Fortier, 05/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function IM_unwrapped=GuidedFloodFill_3D(IM_phase, IM_unwrapped, unwrapped_binary, qualityMap, adjoin, IM_mask,qualityCutoff)


[r_dim, c_dim,z_dim]=size(IM_phase);
IM_mask=logical(IM_mask);

while sum(adjoin(:))~=0             % Loop until there are no more adjoining pixels 
    adjoining_quality=qualityMap + 100.*(1-adjoin);   % Quality values of the adjoining pixels (pad the zero adjoining values with 100 (arbitrary low quality))
    [valueMin,indexMin]=min(adjoining_quality(:));  % Find the highest quality point in the adjoining voxels

    toUnwrap=zeros(size(IM_phase));
    toUnwrap(adjoining_quality<=qualityCutoff*(valueMin))=1;
    toUnwrap=toUnwrap.*adjoin;
    toUnwrap2=find(toUnwrap);
    
    for i=1:length(toUnwrap2)
        [r_adjoin, c_adjoin,z_adjoin]=ind2sub(size(IM_phase),toUnwrap2(i));      % Obtain the coordinates of the first pixel to unwrap
        r_active=r_adjoin(1);
        c_active=c_adjoin(1);
        z_active=z_adjoin(1);
    
        if (IM_mask(r_active, c_active,z_active)==0);  % Unwrap only the region defined by the object mask
           IM_unwrapped(r_active, c_active,z_active)=0;
           adjoin(r_active, c_active,z_active)=0;
        
        % Unwrap the voxel: check every neighbors in an arbitrary order to find a previously unwrapped voxel
        elseif (r_active==r_dim );            
            % search above 
            if unwrapped_binary(r_active-1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active-1, c_active,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search on the right
            elseif c_active~=c_dim &&  unwrapped_binary(r_active, c_active+1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active+1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search on the left
            elseif c_active~=1 && unwrapped_binary(r_active, c_active-1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active-1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search behind
            elseif z_active~=1 && unwrapped_binary(r_active, c_active,z_active-1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active-1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search on the front
            elseif z_active~=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active+1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            else
                adjoin(r_active,c_active,z_active)=0;                                                    % Remove the current active pixel from the adjoin list
            end
        elseif ( r_active==1 );            
        % search below for an adjoining unwrapped phase pixel
            if unwrapped_binary(r_active+1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active+1, c_active,z_active);                                   % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search on the right
            elseif c_active~=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active+1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search on the left
            elseif c_active~=1 && unwrapped_binary(r_active, c_active-1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active-1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search behind
            elseif z_active~=1 && unwrapped_binary(r_active, c_active,z_active-1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active-1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
                
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
           
            % search on the front
            elseif z_active~=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active+1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
              
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
    
            else
                adjoin(r_active,c_active,z_active)=0;                                                    % Remove the current active pixel from the adjoin list
            end
        elseif (c_active==c_dim );            
       % search below for an adjoining unwrapped phase pixel
            if r_active~=r_dim && unwrapped_binary(r_active+1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active+1, c_active,z_active);                                   % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
               
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search above
            elseif r_active~=1 && unwrapped_binary(r_active-1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active-1, c_active,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
          
                % Update the new adjoining pixels:
                if r_active+1>=1 && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search on the left
            elseif c_active~=1 && unwrapped_binary(r_active, c_active-1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active-1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
               
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search behind
            elseif z_active~=1 && unwrapped_binary(r_active, c_active,z_active-1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active-1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
               
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search on the front
            elseif z_active~=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active+1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
     
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            else
                adjoin(r_active,c_active,z_active)=0;                                                    % Remove the current active pixel from the adjoin list
            end
        elseif ( c_active==1);            
        % search below for an adjoining unwrapped phase pixel
            if r_active~=r_dim && unwrapped_binary(r_active+1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active+1, c_active,z_active);                                   % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
               
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search above
            elseif r_active~=1 && unwrapped_binary(r_active-1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active-1, c_active,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
     
                % Update the new adjoining pixels:
                if r_active+1>=1 && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search on the right
            elseif c_active~=c_dim &&  unwrapped_binary(r_active, c_active+1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active+1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
         
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search behind
            elseif z_active~=1 && unwrapped_binary(r_active, c_active,z_active-1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active-1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
         
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search on the front
            elseif z_active~=z_dim &&  unwrapped_binary(r_active, c_active,z_active+1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active+1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
           
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            else
                adjoin(r_active,c_active,z_active)=0;                                                    % Remove the current active pixel from the adjoin list
            end
        elseif (z_active==z_dim );            
        % search below for an adjoining unwrapped phase pixel
            if r_active~=r_dim && unwrapped_binary(r_active+1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active+1, c_active,z_active);                                   % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
        
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            % search above
            elseif r_active~=1 && unwrapped_binary(r_active-1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active-1, c_active,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
          
                % Update the new adjoining pixels:
                if r_active+1>=1 && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            % search on the right
            elseif c_active~=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active+1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
  
                % Update the new adjoining pixels:
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            % search on the left
            elseif c_active~=1 &&  unwrapped_binary(r_active, c_active-1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active-1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
      
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            % search behind
            elseif z_active~=1 && unwrapped_binary(r_active, c_active,z_active-1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active-1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels

                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
            else
                adjoin(r_active,c_active,z_active)=0;                                                    % Remove the current active pixel from the adjoin list
            end
        elseif (z_active==1);            
       % search below for an adjoining unwrapped phase pixel
            if r_active~=r_dim && unwrapped_binary(r_active+1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active+1, c_active,z_active);                                   % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
            
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search above
            elseif r_active~=1 && unwrapped_binary(r_active-1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active-1, c_active,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
   
                % Update the new adjoining pixels:
                if r_active+1>=1 && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search on the right
            elseif c_active~=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active+1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
     
                % Update the new adjoining pixels:
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
            % search on the left
            elseif c_active~=1 && unwrapped_binary(r_active, c_active-1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active-1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels

                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search on the front
            elseif z_active~=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active+1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
       
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
            else
                adjoin(r_active,c_active,z_active)=0;                                                    % Remove the current active pixel from the adjoin list
            end
            
        else
            % search below for an adjoining unwrapped phase pixel
            if unwrapped_binary(r_active+1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active+1, c_active,z_active);                                   % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
      
                % Update the new adjoining pixels:
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1 && adjoin(r_active-1, c_active,z_active)==0
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1 && adjoin(r_active, c_active-1,z_active)==0
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1 && adjoin(r_active-1, c_active,z_active)==0
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1 && adjoin(r_active, c_active,z_active-1)==0
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1 && adjoin(r_active, c_active,z_active+1)==0
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search above
            elseif unwrapped_binary(r_active-1, c_active,z_active)==1
                phase_ref=IM_unwrapped(r_active-1, c_active,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels

                % Update the new adjoining pixels:
                if r_active+1>=1 && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1 && adjoin(r_active+1, c_active,z_active)==0
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1 && adjoin(r_active, c_active-1,z_active)==0
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if c_active+1<=c_dim && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1 && adjoin(r_active, c_active+1,z_active)==0
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1 && adjoin(r_active, c_active,z_active-1)==0
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1 && adjoin(r_active, c_active,z_active+1)==0
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
                
            % search on the right
            elseif unwrapped_binary(r_active, c_active+1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active+1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels

                % Update the new adjoining pixels:
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1 && adjoin(r_active, c_active-1,z_active)==0
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1 && adjoin(r_active-1, c_active,z_active)==0
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1 && adjoin(r_active+1, c_active,z_active)==0
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1 && adjoin(r_active, c_active,z_active-1)==0
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1 && adjoin(r_active, c_active,z_active+1)==0
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
                
            % search on the left
            elseif unwrapped_binary(r_active, c_active-1,z_active)==1
                phase_ref=IM_unwrapped(r_active, c_active-1,z_active);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
         
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1 && adjoin(r_active, c_active+1,z_active)==0
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1 && adjoin(r_active-1, c_active,z_active)==0
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1 && adjoin(r_active+1, c_active,z_active)==0
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1 && adjoin(r_active, c_active,z_active-1)==0
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1 && adjoin(r_active, c_active,z_active+1)==0
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
                
            % search  behind
            elseif unwrapped_binary(r_active, c_active,z_active-1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active-1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
   
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1 && adjoin(r_active, c_active+1,z_active)==0
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1 && adjoin(r_active-1, c_active,z_active)==0
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1 && adjoin(r_active+1, c_active,z_active)==0
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1 && adjoin(r_active, c_active-1,z_active)==0
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active+1<=z_dim && unwrapped_binary(r_active, c_active,z_active+1)==0 && IM_mask(r_active, c_active,z_active+1)==1 && adjoin(r_active, c_active,z_active+1)==0
                    adjoin(r_active, c_active,z_active+1)=1; 
                end
    
            % search on the front
            elseif unwrapped_binary(r_active, c_active,z_active+1)==1
                phase_ref=IM_unwrapped(r_active, c_active,z_active+1);                                       % Obtain the reference unwrapped phase
                p=unwrap([phase_ref IM_phase(r_active, c_active,z_active)]);                             % Unwrap the active pixel
                IM_unwrapped(r_active, c_active,z_active)=p(2);
                unwrapped_binary(r_active, c_active,z_active)=1;                                         % Mark the pixel as unwrapped
                adjoin(r_active, c_active,z_active)=0;                                                   % Remove it from the list of adjoining pixels
        
                % Update the new adjoining pixels:
                if c_active+1>=1 && unwrapped_binary(r_active, c_active+1,z_active)==0 && IM_mask(r_active, c_active+1,z_active)==1 && adjoin(r_active, c_active+1,z_active)==0
                    adjoin(r_active, c_active+1,z_active)=1; 
                end
                if r_active-1>=1 && unwrapped_binary(r_active-1, c_active,z_active)==0 && IM_mask(r_active-1, c_active,z_active)==1 && adjoin(r_active-1, c_active,z_active)==0
                    adjoin(r_active-1, c_active,z_active)=1; 
                end
                if r_active+1<=c_dim && unwrapped_binary(r_active+1, c_active,z_active)==0 && IM_mask(r_active+1, c_active,z_active)==1 && adjoin(r_active+1, c_active,z_active)==0
                    adjoin(r_active+1, c_active,z_active)=1; 
                end
                if c_active-1>=1 && unwrapped_binary(r_active, c_active-1,z_active)==0 && IM_mask(r_active, c_active-1,z_active)==1 && adjoin(r_active, c_active-1,z_active)==0
                    adjoin(r_active, c_active-1,z_active)=1; 
                end
                if z_active-1>=1 && unwrapped_binary(r_active, c_active,z_active-1)==0 && IM_mask(r_active, c_active,z_active-1)==1 && adjoin(r_active, c_active,z_active-1)==0
                    adjoin(r_active, c_active,z_active-1)=1; 
                end
            else
                adjoin(r_active,c_active,z_active)=0;                                                    % Remove the current active pixel from the adjoin list
            end
        end
    end
end
   
