function [data]=read_n4_data_slice_ver5(rawfilename, params, slice, phase, repetition, dataset, echo, acq, partition);
% function [data]=read_n4_data_slice_ver5(rawfilename, params, slice, phase, repetition, dataset, echo, acq, partition);
%
% function to read siemens raw data (Numaris 4)
%      read a single slice, phase, or repetition at a time
%
% inputs:
%     rawfilename  name of raw data file (e.g., 'meas.out')
%     params,  a structure containing:
%          params.protocol_header
%          params.rawheader
%          params.acq.version = 15, 21, 25, or 11
%          params.acq.readoutindicator
%          params.mversion= '6' or '7' (matlab version, defaults to '6')
%     slice
%     phase
%     repetition
%     dataset
%     echo
%     acq
%     partition
%
% output:
%     data(row, col, coil), complex raw data for specified slice, phase, or repetition

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

% check for page algorithm
if ~isfield(params,'mversion')
    params.mversion='6'; % set default for Matlab version 6
end 
if nargin<=5; dataset=1; echo=1; acq=1; end
header=params.rawheader;
readoutindicator=params.acq.readoutindicator;
version=params.acq.version;
if version==15
	kysize=params.protocol_header.sKSpace.lPhaseEncodingLines;
elseif version==21 | version==25
    kysize=params.protocol_header.m_iNoOfFourierLines;
elseif version==11
    kysize=params.protocol_header.iNoOfFourierLines;
end
phase_flag=readoutindicator(:,6)==phase;
slice_flag=readoutindicator(:,3)==slice;
repetition_flag=readoutindicator(:,7)==repetition;
dataset_flag=readoutindicator(:,8)==dataset;
echo_flag=readoutindicator(:,5)==echo;
acq_flag=readoutindicator(:,2)==acq;
partition_flag=readoutindicator(:,4)==partition;
slicereadouts=find((phase_flag+slice_flag+repetition_flag+dataset_flag...
    +echo_flag+acq_flag+partition_flag)==7);

ky=readoutindicator(slicereadouts,1);
readposition=readoutindicator(slicereadouts,12);
ulEvalInfoMask=readoutindicator(slicereadouts,9);
data_size = 2*header.ushSamplesInScan; % real & imag
ncoils=header.ushUsedChannels;
% define ulEvalInfoMask bit masks (from mdh.h include file)
mdh_reflect=base2dec('01000000',16); % define MDH_REFLECT  0x01000000// reflect line

if params.mversion=='7';
    fid=fopen(rawfilename,'r','n');
    tmpdata=zeros(kysize,data_size/2,ncoils,'single');
    for i=1:length(readposition)
        %clear temp
        position=readposition(i);
        fseek(fid,position,-1); %
        for coil=1:ncoils;
            temp(:,coil) = fread(fid,data_size,'single => single');
            fseek(fid,128,0); % skip mdh header
        end
        temp2 =complex(temp(1:2:data_size,:), temp(2:2:data_size,:)); % separate the real and imaginary component
        reflect_flag=bitand(ulEvalInfoMask(i),mdh_reflect)>0;
        if reflect_flag==1; temp2=temp2(end:-1:1,:);end
        tmpdata(ky(i),:,:)=temp2;
    end
    fclose(fid);
    % zero values of cut-off readout
    data=zeros(size(tmpdata),'single');
    kxmin=header.sCutOffData.ushPre + 1; % /* write ushPre zeros at line start */
    kxmax=size(data,2) - header.sCutOffData.ushPost; % /* write ushPost zeros at line end */
    data(:,kxmin:kxmax,:)=tmpdata(:,kxmin:kxmax,:);
else
    fid=fopen(rawfilename,'r','n');
    tmpdata=zeros(kysize,data_size/2,ncoils);
    for i=1:length(readposition)
        %clear temp
        position=readposition(i);
        fseek(fid,position,-1); %
        for coil=1:ncoils;
            temp(:,coil) = fread(fid,data_size,'single');
            fseek(fid,128,0); % skip mdh header
        end
        temp2 =complex(temp(1:2:data_size,:), temp(2:2:data_size,:)); % separate the real and imaginary component
        reflect_flag=bitand(ulEvalInfoMask(i),mdh_reflect)>0;
        if reflect_flag==1; temp2=temp2(end:-1:1,:);end
        tmpdata(ky(i),:,:)=temp2;
    end
    fclose(fid);
    % zero values of cut-off readout
    data=zeros(size(tmpdata));
    kxmin=header.sCutOffData.ushPre + 1; % /* write ushPre zeros at line start */
    kxmax=size(data,2) - header.sCutOffData.ushPost; % /* write ushPost zeros at line end */
    data(:,kxmin:kxmax,:)=tmpdata(:,kxmin:kxmax,:);
end    
    
return
