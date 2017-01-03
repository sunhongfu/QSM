function [header, readoutindicator, readoutindicator_full] = read_n4_indicators_ver3(rawfilename, version)
% function [header, readoutindicator, readoutindicator_full] = read_n4_indicators_ver2(rawfilename, version);
%
% function to read siemens raw data file readout indicator (Numaris 4)
%
% input:
%     rawfilename  name of raw data file (e.g., 'meas.out')
%     version  15, 21, or 25 (default=21)
% output:
%     readoutindicator is a matrix with a row vector for each readout containing:
%              [kyline,acquisition,slice,partition,echo,phase,repetition,dataset,...
%                          ulEvalInfoMask,ulTimeStamp,ulPMUTimeStamp, readoutposition]
%              for example: plot(readoutindicator(:,1),'.') plots the phase encode order
%                      or:  plot(readoutindicator(:,10),'.') plots the time
%                                                           stamp of each
%                                                           readout (2.5 ms ticks)
%     readoutindicator_full includes all readouts
%     readoutindicator excludes prescan noise readouts
%     header is the 1-st mdh header (for 1-st readout following pre-scan)

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

% --- 2013-10-31 Kelvin: Version 3 changes -------------------------------------
% DMALength may not accurately reflect the number of bytes of each raw data
% readout, and therefore total_readouts may not be accurate. Instead, simply
% continue reading raw data lines until no more data remains

if nargin==1; version=21; end % set default version

fid=fopen(rawfilename,'r','n');
fseek(fid,32,-1); % skip over the 32 byte pre-header
header1 = read_n4_header(fid,version); % read first header to get DMALength (i.e. the number of bytes of raw
                            % data per readout plus the 128 byte header
fseek(fid,0,1);
total_length = ftell(fid); % seek to the end of the data file to find the total length
data_size_per_readout = double(header1.ulDMALength); % (samples in scan * 8 + mdh size) * ncoils
total_readouts = floor((total_length - 32)/data_size_per_readout);  % total number of readouts is calculated
fseek(fid,32,-1); % Reposition the fid index to just skip over the 32 byte pre-header
ncoils=header1.ushUsedChannels;
ushSamplesInScan=header1.ushSamplesInScan;

readoutindicator_full=zeros(total_readouts,12);
if version==15; n1=0; elseif version==21 | version==25 | version == 11; n1=4; end
n2=128-9*2-n1-7*4;

% for j = 1:total_readouts;  % Kelvin
j = 1;
while(1)
    headerdata=fread(fid,7,'uint32');
    fseek(fid,n1,0);
    counters=1+fread(fid,9,'uint16');
		
		% Kelvin
		if isempty(counters)
			break
		end

    fseek(fid,n2,0);
    position=ftell(fid);% file position for beginning of data for j-th readout
    offset=128*(ncoils-1)+ushSamplesInScan*4*2*ncoils;
    fseek(fid,offset,0);
    readoutindicator_full(j,:)=[...
        counters(1),...
		counters(2),...
		counters(3),...
		counters(4),...
		counters(5),...
		counters(6),...
		counters(7),...
		counters(8),...
        headerdata(6),...
        headerdata(4),...
        headerdata(5),...
        position];
			
		j = j+1;  % Kelvin
end % readout loop

readoutindicator=check_valid_data(readoutindicator_full);

% position of 1-st valid mdh header (after pre-scan noise)
mdh_header1_position=readoutindicator(1,12)-128;
fseek(fid,mdh_header1_position,-1);
[header] = read_n4_header(fid, version);


fclose (fid);
return

function [header] = read_n4_header(fid, version);
% function [header] = read_n4_header(fid,version);
%
% function to read siemens raw data file header (version Numaris 4)
%
% input:
%     fid     file handle, e.g.   fid=fopen(rawfilename,'r','n');
%     version 15, 21, or 25
% output:
%     header  structure containing MDH header fields.

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

if nargin==1; version=21; end % set default

if version ==15 | version ==21
	header.ulDMALength         = fread(fid,1,'uint32'); % DMA length [bytes]           
elseif version==25 | version==11 
    header.ulDMALength         = fread(fid,1,'uint32'); % bit  0..27: DMA length [bytes]
														% bit 28..31: pci_rx enable flags
end
header.lMeasUID            = fread(fid,1,'int32'); % measurement user ID
header.ulScanCounter       = fread(fid,1,'uint32');  % scan counter [1...]
header.ulTimeStamp         = fread(fid,1,'uint32');  % time stamp [2.5 ms ticks since 00:00]
header.ulPMUTimeStamp      = fread(fid,1,'uint32'); % PMU time stamp [2.5 ms ticks since last trigger]
if version==15
	header.ulEvalInfoMask      = fread(fid,1,'uint32');  % evaluation info mask
elseif version==21 | version==25
    header.ulEvalInfoMask      = fread(fid,2,'uint32');  % evaluation info mask
end
header.ushSamplesInScan    = fread(fid,1,'uint16');  % # of samples acquired in scan
header.ushUsedChannels     = fread(fid,1,'uint16'); % # of channels used in scan

% loop counters
header.sLC.ushLine         = 1+fread(fid,1,'uint16'); % line index                   
header.sLC.ushAcquisition  = 1+fread(fid,1,'uint16'); % acquisition index            
header.sLC.ushSlice        = 1+fread(fid,1,'uint16'); % slice index                  
header.sLC.ushPartition    = 1+fread(fid,1,'uint16'); % partition index              
header.sLC.ushEcho         = 1+fread(fid,1,'uint16'); % echo index                   
header.sLC.ushPhase        = 1+fread(fid,1,'uint16'); % phase index                  
header.sLC.ushRepetition   = 1+fread(fid,1,'uint16'); % measurement repeat index     
header.sLC.ushSet          = 1+fread(fid,1,'uint16'); % set index                    
header.sLC.ushSeg          = 1+fread(fid,1,'uint16'); % segment index  (for TSE)     
if version==15
	header.sLC.ushFree         = fread(fid,1,'ushort'); %  free loop counter            
elseif version==21 | version==25;
    header.sLC.ushIda =fread(fid,1,'ushort');  % IceDimension a index        
    header.sLC.ushIdb =fread(fid,1,'ushort');  % IceDimension b index        
    header.sLC.ushIdc =fread(fid,1,'ushort');  % IceDimension c index        
    header.sLC.ushIdd =fread(fid,1,'ushort');  % IceDimension d index        
    header.sLC.ushIde =fread(fid,1,'ushort');  % IceDimension e index        
end
% cut-off values  
header.sCutOffData.ushPre  = fread(fid,1,'uint16'); % write ushPre zeros at line start 
header.sCutOffData.ushPost = fread(fid,1,'uint16'); % write ushPost zeros at line end  
header.ushKSpaceCentreColumn   = fread(fid,1,'uint16');  % centre of echo
if version ==15 | version==21
	header.ushDummy            = fread(fid,1,'uint16');  % for swapping
elseif version ==25
    header.ushCoilSelect       = fread(fid,1,'uint16');  % Bit 0..3: CoilSelect
end
header.fReadOutOffcentre       = fread(fid,1,'single');  % ReadOut offcenter value
header.ulTimeSinceLastRF       = fread(fid,1,'int32');   % Sequence time stamp since last RF pulse           
header.ushKSpaceCentreLineNo   = fread(fid,1,'uint16');  % number of K-space centre line
header.ushKSpaceCentrePartitionNo = fread(fid,1,'uint16');  % number of K-space centre partition
if version==15
    fseek(fid,28,0); % skip over free data buffer
elseif version==21 | version==25
    header.aushIceProgramPara=fread(fid,4,'uint16');
    header.aushFreePara=fread(fid,4,'uint16');
end
% slice data (28 bytes)
header.sSD.sVector.flSag =  fread(fid,1,'single');%  slice position vector        
header.sSD.sVector.flCor =  fread(fid,1,'single');%  slice position vector        
header.sSD.sVector.flTra =  fread(fid,1,'single');%  slice position vector        
header.sSD.aflQuaternion.a =  fread(fid,1,'single'); % rotation matrix as quaternion
header.sSD.aflQuaternion.b =  fread(fid,1,'single'); % rotation matrix as quaternion
header.sSD.aflQuaternion.c =  fread(fid,1,'single'); % rotation matrix as quaternion
header.sSD.aflQuaternion.d =  fread(fid,1,'single'); % rotation matrix as quaternion
header.ulChannelId = fread(fid,1,'int32');  % channel Id must be the last parameter

return

function readoutindicator_valid_data=check_valid_data(readoutindicator);
% function readoutindicator_valid_data=check_valid_data(readoutindicator);
%
% function to delete rows (readouts) from readoutindicator which correspond
% to noiseadjscan or reference scan readouts

% define ulEvalInfoMask bit masks (from mdh.h include file)
mdh_noiseadjscan=base2dec('02000000',16);% define MDH_NOISEADJSCAN  0x02000000 // noise adjust scan
mdh_extrarefscan=base2dec('00002000',16);% define MDH_PPAEXTRAREFSCAN  0x00002000 // additional scan for PPA reference line/partition
mdh_phasecorrection=base2dec('00200000',16); % define MDH_PHASCOR 0x00200000 // phase correction data    


noiseadjscan_flag=bitand(readoutindicator(:,9),mdh_noiseadjscan)>0;
extrarefscan_flag=bitand(readoutindicator(:,9),mdh_extrarefscan)>0;
phasecorrection_flag=bitand(readoutindicator(:,9),mdh_phasecorrection)>0;
valid_flag=(noiseadjscan_flag==0) & (extrarefscan_flag==0) &(phasecorrection_flag==0);
validreadouts=find(valid_flag==1);
readoutindicator_valid_data=readoutindicator(validreadouts,:);

return

