function IM3D=reshapemat(cpxdata,Pars,shiftbyte,ebytes)

% reshapes raw scanner output into 3D matrix

if ebytes==4      
    shiftbyte1=-1;
else
    shiftbyte1=0;
end

if Pars.nv2>1 % If data is 3D, follow method to remove unwanted data points
    
    % In 3D raw data, found there is an extra 7 points at the start of
    % every slice that need to be disposed of
    
    if Pars.RCVRS_==4  %4 RCVRS
        % 3D 4 RCVR data is ordered differently than 2D, so requires
        % different matrix reshaping. 
        cpxdata=reshape(cpxdata,[],Pars.nv2*Pars.RCVRS_);
        % cpxdata=cpxdata(8:end,:);
        cpxdata=cpxdata(15:end,:);  %%% HONGFU added

        %In 3D data, 4 RCVRS 
        cpxdata1=cpxdata(:,1:4:end); % RCVR 1
        cpxdata2=cpxdata(:,2:4:end); % RCVR 2
        cpxdata3=cpxdata(:,3:4:end); % RCVR 3
        cpxdata4=cpxdata(:,4:4:end); % RCVR 4
                               
%         %First RCVR
        cpxdata1=reshape(cpxdata1,1,[]);
        cpxdata1=cpxdata1';
        cpxdata1=reshape(cpxdata1,2,[]);
        cpxdata1=rot90(cpxdata1);
        cpxdata1=complex(cpxdata1(:,1),cpxdata1(:,2));
        IM3D.RCVR1=reshape(cpxdata1,Pars.np/2,Pars.nv,[]);
        clear cpxdata1
        IM3D.RCVR1=single(IM3D.RCVR1);        
         
%         %2nd RCVR
        cpxdata2=reshape(cpxdata2,1,[]);
        cpxdata2=cpxdata2';
        cpxdata2=reshape(cpxdata2,2,[]);
        cpxdata2=rot90(cpxdata2);
        cpxdata2=complex(cpxdata2(:,1),cpxdata2(:,2));
        IM3D.RCVR2=reshape(cpxdata2,Pars.np/2,Pars.nv,[]);
        clear cpxdata2
        IM3D.RCVR2=single(IM3D.RCVR2);  
        
%         %3rd RCVR
        cpxdata3=reshape(cpxdata3,1,[]);
        cpxdata3=cpxdata3';
        cpxdata3=reshape(cpxdata3,2,[]);
        cpxdata3=rot90(cpxdata3);
        cpxdata3=complex(cpxdata3(:,1),cpxdata3(:,2));
        IM3D.RCVR3=reshape(cpxdata3,Pars.np/2,Pars.nv,[]);
        clear cpxdata3
        IM3D.RCVR3=single(IM3D.RCVR3);    
    
%         %4th RCVR
        cpxdata4=reshape(cpxdata4,1,[]);
        cpxdata4=cpxdata4';
        cpxdata4=reshape(cpxdata4,2,[]);
        cpxdata4=rot90(cpxdata4);
        cpxdata4=complex(cpxdata4(:,1),cpxdata4(:,2));
        IM3D.RCVR4=reshape(cpxdata4,Pars.np/2,Pars.nv,[]);
        clear cpxdata4
        IM3D.RCVR4=single(IM3D.RCVR4);
  
    else  %1 RCVR
        cpxdata1=reshape(cpxdata,[],Pars.nv2); % this could apply to the 2D data without causing any harm
        clear cpxdata
        cpxdata1=cpxdata1(8:end,:);
        cpxdata1=reshape(cpxdata1,1,[]);
        cpxdata1=cpxdata1';
        cpxdata1=reshape(cpxdata1,2,[]);
        cpxdata1=rot90(cpxdata1);
        cpxdata1=complex(cpxdata1(:,1),cpxdata1(:,2));
        IM3D.RCVR1=reshape(cpxdata1,Pars.np/2,Pars.nv,[]);
        clear cpxdata1
        IM3D.RCVR1=single(IM3D.RCVR1);

    end
    
else  %If data is 2D
    
    if strcmp(Pars.seqfil,'ge3d')  % Amir's 2D sequence
        if ebytes==4  % Should this only be for 1 RCVR???
            % cpxdata=cpxdata(1:end-7,:);
            cpxdata=cpxdata(8:end,:);  %%% HONGFU'S fix
        end

        if Pars.RCVRS_==1
            cpxdata=single(cpxdata);
            cpxdata=circshift(cpxdata,shiftbyte1+shiftbyte);
        end  
    end
    
    if Pars.RCVRS_==4   % 4 receivers
        % 4 in the next line = 4 receivers
        cpxdata=reshape(cpxdata,[],4);
        cpxdata=cpxdata(1:end-14,:);
        cpxdata1=cpxdata(:,1); cpxdata2=cpxdata(:,2); cpxdata3=cpxdata(:,3); cpxdata4=cpxdata(:,4); 
        clear cpxdata
    
        %Reshaping to real and imaginary (complex)
        cpxdata1=reshape(cpxdata1,2,[]); %%%%This line is causing problems for 1 RCVR!!!
        cpxdata1=single(cpxdata1);
        cpxdata1=rot90(cpxdata1);  %GIVING MEMORY ERRORS!!!
        cpxdata1=complex(cpxdata1(:,1),cpxdata1(:,2));
        
        cpxdata2=reshape(cpxdata2,2,[]); %%%%This line is causing problems for 1 RCVR!!!
        cpxdata2=single(cpxdata2);
        cpxdata2=rot90(cpxdata2);  %GIVING MEMORY ERRORS!!!
        cpxdata2=complex(cpxdata2(:,1),cpxdata2(:,2));
        
        cpxdata3=reshape(cpxdata3,2,[]); %%%%This line is causing problems for 1 RCVR!!!
        cpxdata3=single(cpxdata3);
        cpxdata3=rot90(cpxdata3);  %GIVING MEMORY ERRORS!!!
        cpxdata3=complex(cpxdata3(:,1),cpxdata3(:,2));
        
        cpxdata4=reshape(cpxdata4,2,[]); %%%%This line is causing problems for 1 RCVR!!!
        cpxdata4=single(cpxdata4);
        cpxdata4=rot90(cpxdata4);  %GIVING MEMORY ERRORS!!!
        cpxdata4=complex(cpxdata4(:,1),cpxdata4(:,2));
        
%         if Pars.nv2>1 % Data is 3D
%             IM3D.RCVR1=reshape(f1,Pars.np/2,Pars.nv,Pars.nv2); clear f1
%             IM3D.RCVR2=reshape(f2,Pars.np/2,Pars.nv,Pars.nv2); clear f2
%             IM3D.RCVR3=reshape(f3,Pars.np/2,Pars.nv,Pars.nv2); clear f3
%             IM3D.RCVR4=reshape(f4,Pars.np/2,Pars.nv,Pars.nv2); clear f4    
%         else % Data is 2D
            IM3D.RCVR1=reshape(cpxdata1,Pars.np/2,Pars.NS,Pars.nv); clear cpxdata1
            IM3D.RCVR2=reshape(cpxdata2,Pars.np/2,Pars.NS,Pars.nv); clear cpxdata2
            IM3D.RCVR3=reshape(cpxdata3,Pars.np/2,Pars.NS,Pars.nv); clear cpxdata3
            IM3D.RCVR4=reshape(cpxdata4,Pars.np/2,Pars.NS,Pars.nv); clear cpxdata4
%         end
    
        IM3D.RCVR1=single(IM3D.RCVR1); 
        IM3D.RCVR2=single(IM3D.RCVR2); 
        IM3D.RCVR3=single(IM3D.RCVR3); 
        IM3D.RCVR4=single(IM3D.RCVR4); 

    elseif Pars.RCVRS_==1  % 1 receiver
        
        if strcmp(Pars.seqfil,'gems') % gems/varian 2D sequence as opposed to Amir's 2D sequence
            cpxdata=cpxdata(8:end);
        end
        
        % Reshaping to real and imaginary (complex)
        cpxdata=reshape(cpxdata,2,[]); %%%%This line is causing problems for 1 RCVR!!!
        % for 3D 1 RCVR data too
        cpxdata=single(cpxdata);
        cpxdata=rot90(cpxdata);  %GIVING MEMORY ERRORS!!!
        cpxdata=complex(cpxdata(:,1),cpxdata(:,2));
        
        % cpxdata=cpxdata(1:end-7,:);   %%% HONGFU commented out
        f1=cpxdata(:,1); 

        clear cpxdata
%         if Pars.nv2>1 %3D data
%             IM3D.RCVR1=reshape(f1,Pars.np/2,Pars.nv,[]); clear f1
%         else   % 2D data
            IM3D.RCVR1=reshape(f1,Pars.np/2,Pars.NS,Pars.nv); clear f1
%         end
        IM3D.RCVR1=single(IM3D.RCVR1);
    end
    
end