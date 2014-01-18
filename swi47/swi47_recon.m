function [img,Pars] = swi47_recon(path_in)

if ~ exist('path_in','var') || isempty(path_in)
    path_in = [pwd '/ge3d__01.fid'];
end

file=[path_in,'/fid'];
parfil=[path_in,'/procpar'];
Pars = Par_Read_Struct(parfil);
% Generalize to both Amir and Varian 2D acquisition
if Pars.ns < 2
Pars.np=Pars.np/Pars.NS;
end
Pars.ns=Pars.NS;


% READING THE VARIAN FILE (RAW DATA)
[fid,errmsg] = fopen(file,'r','ieee-be');
if fid < 0
ErrMsg = [file 13 10 errmsg];
error(ErrMsg);
end
[datafilehead, count] = fread(fid,8,'int32');
% struct datafilehead
nblocks = datafilehead(1);	%number of blocks in file
ntraces = datafilehead(2);	% number of traces per block
Pars.npX = datafilehead(3);	% number of elements per trace (Real + Imag)
ebytes = datafilehead(4);	% number of bytes per element
tbytes = datafilehead(5);	% number of bytes per trace
bbytes = datafilehead(6);	% number of bytes per block
vers_id = rem(datafilehead(7),32768);	% software version and file_id status bits
status = floor(datafilehead(7)/32768);	% status of whole file
nbdatafileheads = datafilehead(8);	% number of block datafileheads
filehead = [nblocks, ntraces, Pars.np, ebytes, tbytes, bbytes, vers_id, status, nbdatafileheads];
% struct datablockhead
[scale, count] = fread(fid,1,'int16');	 % scaling factor
[status, count] = fread(fid,1,'int16'); % status of data in block
[index, count] = fread(fid,1,'int16');  % block index
[mode, count] = fread(fid,1,'int16');   % mode of data in block
[ctcount, count] = fread(fid,1,'int32');   % ct value for FID
[lpval, count] = fread(fid,1,'float32');   % F2 left phase in phasefile
[rpval, count] = fread(fid,1,'float32');   % F2 right phase in phasefile
[lvl, count] = fread(fid,1,'float32');     % F2 level drift correction
[tlt, count] = fread(fid,1,'float32');     % F2 tilt drift correction
blockhead = [scale, status, index, mode, ctcount,lpval,rpval,lvl,tlt];
stat= fseek(fid,32, 'bof');	% Position file pointer after datafilehead
datasize = [700+ntraces*Pars.np*nblocks];
% Detect the variable size and read the file in
if ebytes==2
[cpxdata, count] = fread(fid,inf ,'*int16');
elseif ebytes==4
[cpxdata, count] = fread(fid,inf ,'*int32');
else
k=strcat('ebytes = ', num2str(ebytes), ' this element size is different from the ')
k='two acceptable types .. revise me'
R = iPars.nput('Stop program now to decide what to do (press CTRL C) ')
end
fclose(fid);


% Reshape raw scanner output to 3D matrix
IM3D=reshapemat(cpxdata,Pars,1,ebytes);
IM3D_RCVR1=IM3D.RCVR1;
if Pars.RCVRS_==4  % 4 RCVRS
IM3D_RCVR2=IM3D.RCVR2;
IM3D_RCVR3=IM3D.RCVR3;
IM3D_RCVR4=IM3D.RCVR4;
end
clear cpxdata IM3D 


% REMOVE DC SHIFTS/ARTIFACTS
IM3D_RCVR1=unDCshift(IM3D_RCVR1,Pars);
if Pars.RCVRS_==4
IM3D_RCVR2=unDCshift(IM3D_RCVR2,Pars);
IM3D_RCVR3=unDCshift(IM3D_RCVR3,Pars);
IM3D_RCVR4=unDCshift(IM3D_RCVR4,Pars);
end


% RECONSTRUCT PARTIAL READOUT DATA
if Pars.RO_CVRD<1   % If fraction of readout is less than 1 (partial readout)
IM3D_RCVR1=partial_RO(IM3D_RCVR1,Pars);
if Pars.RCVRS_==4
IM3D_RCVR2=partial_RO(IM3D_RCVR2,Pars);
IM3D_RCVR3=partial_RO(IM3D_RCVR3,Pars);
IM3D_RCVR4=partial_RO(IM3D_RCVR4,Pars);
end
end


% HANNING FILTER AND WEIGHT RCVRS
Hanwindow3D=Hanning_win(Pars);
if Pars.RCVRS_==1
IM3D_RCVR1=IM3D_RCVR1.*Hanwindow3D;
clear Hanwindow3D
elseif Pars.RCVRS_==4
IM3D_RCVR1=IM3D_RCVR1.*Hanwindow3D;
IM3D_RCVR2=IM3D_RCVR2.*Hanwindow3D;
IM3D_RCVR3=IM3D_RCVR3.*Hanwindow3D;
IM3D_RCVR4=IM3D_RCVR4.*Hanwindow3D;
clear Hanwindow3D
% multiply each receiver's output by its weighting factor
IM3D_RCVR1=IM3D_RCVR1*1; %was 1
IM3D_RCVR2=IM3D_RCVR2*0.97; % was 0.97
IM3D_RCVR3=IM3D_RCVR3*1.01; % was 1.01;     4/2.4
IM3D_RCVR4=IM3D_RCVR4*1.12; %was 1.12;       4/2.4
end


% If data is 3D, add a ramp to shift k-space to proper location
% Also permute matrix so dimensions are in the same order as 2D (for
% purpose of filtering/viewing) i.e. (x,z,y)
IM3D_RCVR1=kramp(Pars,IM3D_RCVR1);
if Pars.RCVRS_==4
IM3D_RCVR2=kramp(Pars,IM3D_RCVR2);
IM3D_RCVR3=kramp(Pars,IM3D_RCVR3);
IM3D_RCVR4=kramp(Pars,IM3D_RCVR4);
end


% FOURIER TRANSFORM TO IMAGE SPACE
% Slice (Pars.NS ) direction in this case is the second dimePars.nsion of the matrix
if Pars.RCVRS_==4
IM3D_RCVR1=fft(single(fftshift(IM3D_RCVR1)),[],1);
IM3D_RCVR1=fft(IM3D_RCVR1,[],3);
IM3D_RCVR1=fftshift(IM3D_RCVR1);
IM3D_RCVR2=fft(single(fftshift(IM3D_RCVR2)),[],1);
IM3D_RCVR2=fft(IM3D_RCVR2,[],3);
IM3D_RCVR2=fftshift(IM3D_RCVR2);
IM3D_RCVR3=fft(single(fftshift(IM3D_RCVR3)),[],1);
IM3D_RCVR3=fft(IM3D_RCVR3,[],3);
IM3D_RCVR3=fftshift(IM3D_RCVR3);
IM3D_RCVR4=fft(single(fftshift(IM3D_RCVR4)),[],1);
IM3D_RCVR4=fft(IM3D_RCVR4,[],3);
IM3D_RCVR4=fftshift(IM3D_RCVR4);
if Pars.nv2>1  % 3D data - also needs to be FT'd in the 3rd dimension
IM3D_RCVR1=fftshift(IM3D_RCVR1);
IM3D_RCVR1=fft(IM3D_RCVR1,[],2);
IM3D_RCVR1=fftshift(IM3D_RCVR1);
IM3D_RCVR2=fftshift(IM3D_RCVR2);
IM3D_RCVR2=fft(IM3D_RCVR2,[],2);
IM3D_RCVR2=fftshift(IM3D_RCVR2);
IM3D_RCVR3=fftshift(IM3D_RCVR3);
IM3D_RCVR3=fft(IM3D_RCVR3,[],2);
IM3D_RCVR3=fftshift(IM3D_RCVR3);
IM3D_RCVR4=fftshift(IM3D_RCVR4);
IM3D_RCVR4=fft(IM3D_RCVR4,[],2);
IM3D_RCVR4=fftshift(IM3D_RCVR4);
end
elseif Pars.RCVRS_==1  % 1 receiver
% Slice (Pars.NS ) direction in this case is the second dimension of the matrix
IM3D_RCVR1=fft(single(fftshift(IM3D_RCVR1)),[],1);
IM3D_RCVR1=fft(IM3D_RCVR1,[],3);
IM3D_RCVR1=fftshift(IM3D_RCVR1);
if Pars.nv2>1  % 3D data - also needs to be FT'd in the 3rd dimension
IM3D_RCVR1=fftshift(IM3D_RCVR1);
IM3D_RCVR1=fft(IM3D_RCVR1,[],2);
IM3D_RCVR1=fftshift(IM3D_RCVR1);
end
end


% PUT MATRIX DIMENSIONS INTO PROPER ORDER (X,Y,Z) (Most 2D data has
% dimensions ordered differently - x,z,y)
IM3D_RCVR1=PermuteXYZ(IM3D_RCVR1,Pars);
if Pars.RCVRS_==4
IM3D_RCVR2=PermuteXYZ(IM3D_RCVR2,Pars);
IM3D_RCVR3=PermuteXYZ(IM3D_RCVR3,Pars);
IM3D_RCVR4=PermuteXYZ(IM3D_RCVR4,Pars);
end

img = cat(4,IM3D_RCVR1,IM3D_RCVR2,IM3D_RCVR3,IM3D_RCVR4);