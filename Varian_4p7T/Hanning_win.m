function Hanwindow3D=Hanning_win(Pars)

% Produce Hanning window and apply to raw data

% Produced by Sandra Meyers/Amir Eissa
% June 19/09

% Input: Parameters variable


% Hanning 25% on each side of each data vector (Room to change the Hanning
% percentage .... also need to add the 3rd Hanning vector in the 3D imaging
% case
HanWin=hann(round((Pars.np/2)/2)); % x Hanning
HanWin2=hann(round(Pars.nv/2)); % y Hanning

HanWin=[HanWin(1:round(max(size(HanWin))/2))' ones(1,Pars.np/2-round(max(size(HanWin))/2)*2) HanWin(round(max(size(HanWin))/2)+1:end)'];
HanWin2=[HanWin2(1:round(max(size(HanWin2))/2))' ones(1,Pars.nv-round(max(size(HanWin2))/2)*2) HanWin2(round(max(size(HanWin2))/2)+1:end)'];

%if Pars.nv2>1   % If data is 3D, add a 3rd Hann window in z direction
    %HanWin3=hann(round(Pars.nv2/2));
    %HanWin3=[HanWin3(1:round(max(size(HanWin3))/2))' ones(1,Pars.nv2-round(max(size(HanWin3))/2)*2) HanWin3(round(max(size(HanWin3))/2)+1:end)'];
    %[Han1,Han2,Han3]=meshgrid(HanWin2,HanWin1,HanWin3);
    % Note that the HanWins have to be this order so that meshgrid will order
    % Hans properly as x,y,z
    %Han4=Han1.*Han2.*Han3;
%else     % 2D data
    [Han1,Han2]=meshgrid(HanWin,HanWin2);   
    Han3=Han1.*Han2;  % Check matrix orientation (M3 should have same size and orientation as IM3D_RCVR1
    Hanwindow3D=repmat(Han3,[1 1 max(Pars.NS,Pars.nv2)]); % Note that using this max generalizes to 2D OR 3D data
    if Pars.nv2>1
        Hanwindow3D=permute(Hanwindow3D,[2 1 3]);  %If data is 3D, dimensions fall in a different order
    else
        Hanwindow3D=shiftdim(Hanwindow3D,1);
    end
%end

