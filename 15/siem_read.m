function [mrdata,read_status] = siem_read(pathfile)

% siem_read.m
%
% function to read siemens raw data files

%     *************************************************
%     *  Richard Thompson  (thompsor@nhlbi.nih.gov)   *
%     *  Laboratory for Cardiac Energetics            *
%     *  NIH NHLBI                                    *
%     *************************************************  
%
% June 2013 Alan modified for 3D imaging too

if nargin < 1
  [FILENAME, PATHNAME] = uigetfile('*.out', 'Please select the file to open');
  read_status = 1;
  if isequal(FILENAME,0)|isequal(PATHNAME,0) 
    read_status = 0;
    mrdata = [];
    return;
  end
  
  FILENAMEfull=[PATHNAME FILENAME]; % header name and path
else
  PATHNAME = cell2mat(pathfile(1));
  FILENAME = cell2mat(pathfile(2));
  FILENAMEfull=strcat(PATHNAME,FILENAME); % header name and path
end



fid = fopen(FILENAMEfull,'r'); % open the raw data file

if fid > 0
    
    ll = length(FILENAMEfull);
    filename = FILENAMEfull;
    tempname = FILENAME(1:end-8);        
    tempname2 = strcat(tempname,'meas.asc');
    %tempname2 = strcat(tempname,'MrProt.asc');
    rawfilename=FILENAMEfull;
    asciifilename =[PATHNAME tempname2];
    
    % read raw file
    params.acq.path = PATHNAME;
    params.acq.rawfilename=rawfilename;
    params.acq.asciifilename = asciifilename;
    
    params.protocol_header.ulVersion = read_protocol_header_fieldname(asciifilename,'ulVersion');
    
    switch params.protocol_header.ulVersion
        case ' 0x3e9'
            params.acq.version=15;
        case ' 0x7da'
            params.acq.version=21;
        case ' 0xbee332'
            params.acq.version=25;
        case ' 0x1421cf5'
            params.acq.version=11;
    end
    
    params.acq.slice = read_protocol_header_fieldname(asciifilename,'sSliceArray.lSize');
    
    % a fix to arbitraily set the version to 21 - not sure if this is
    % appropriate
    if ~isfield(params.acq,'version')
      params.acq.version = 21;
    end
    
    [rawheader, params.acq.readoutindicator, params.acq.readoutindicator_full] = ...
        read_n4_indicators_ver2(rawfilename, params.acq.version);
    readoutindicator=params.acq.readoutindicator;
    
    params.protocol_header.sKSpace.lPhaseEncodingLines = read_protocol_header_fieldname(asciifilename,'sKSpace.lPhaseEncodingLines');
    params.acq.rawfilename=rawfilename;
    protocol_header = read_protocol_header(asciifilename);
    params.protocol_header=protocol_header;
    % params.protocol_header.sKSpace.lPhaseEncodingLines = 1; % Hongfu comment out
    params.rawheader=rawheader;
    params.acq.slice=1; phase=1; repetition=1; dataset = 1; echo1=1; acq = 1;
    [data]=read_n4_data_slice_ver5(rawfilename, params, params.acq.slice, phase, repetition, dataset, echo1, acq, 1);
    params.acq.kxsize=size(data,2);
    params.acq.kxcenter = params.acq.kxsize/2;
    params.window.param_ref = 0;
    params.acq.xsize = params.acq.kxsize;
    params.acq.ysize = params.acq.kxsize;
    params.window.param_data = 0;
    params.homodyne.flag = 1;
    params.homodyne.Niterations = 2;
    params.acq.kxsize2=read_protocol_header_fieldname(asciifilename,'sKSpace.lBaseResolution');
    params.acq.kysize=size(data,1);% protocol_header.sKSpace.lPhaseEncodingLines;
    params.acq.ncoils=params.rawheader.ushUsedChannels;
    params.acq.nacquisitions=max(params.acq.readoutindicator(:,2));
    params.acq.nslices=max(params.acq.readoutindicator(:,3));
    params.acq.npartitions=max(params.acq.readoutindicator(:,4));
    params.acq.nechos=max(params.acq.readoutindicator(:,5));
    params.acq.nphases=max(params.acq.readoutindicator(:,6));
    params.acq.nsets=max(params.acq.readoutindicator(:,8));
    params.acq.nrepetitions=max(params.acq.readoutindicator(:,7));
    
    if params.acq.kxsize/params.acq.kxsize2==2
        params.acq.readoutoversampling=1;
    else
        params.acq.readoutoversampling=0;
    end
    
    for phase=1:params.acq.nphases
        tmp1=readoutindicator(:,6)==phase;
        tmp2=readoutindicator(:,3)==1; % slice;
        tmp3=readoutindicator(:,7)==1; %repetition;
        tmp4=readoutindicator(:,8)==1; %dataset;
        slicereadouts=find((tmp1+tmp2+tmp3+tmp4)==4);
        ky{phase}=sort(readoutindicator(slicereadouts,1));
%        ky{phase}=sort(readoutindicator(readoutindicator(:,6)==phase,1))';
    end
    
    params.view=ky;
    
    %window=repmat(gaussian_window([params.acq.kxsize params.acq.kxsize], 1),[1 1 params.acq.ncoils]);
    %params.acq.nphases = 4;
    
    params.acq.kysize;
    mrdata.raw=zeros(params.acq.kysize,params.acq.kxsize,params.acq.npartitions,params.acq.ncoils,...
    params.acq.nphases,params.acq.nsets,params.acq.nacquisitions,params.acq.nechos,...
    params.acq.nrepetitions,params.acq.nslices,'single');
    
%    h = waitbar(0,'Please wait...');
    
    total_readouts = params.acq.nphases*params.acq.nsets*params.acq.nechos*params.acq.nacquisitions*params.acq.nrepetitions;
    counter = 0;
    for cardiac_phase=1:params.acq.nphases
%        disp(['processing cardiac phase: ',num2str(cardiac_phase),' of ',num2str(params.acq.nphases),' phases'])
        tic
        for slice = 1:params.acq.nslices
          for repetition = 1:params.acq.nrepetitions
            for dataset=1:params.acq.nsets
              for echo1=1:params.acq.nechos
                for acq=1:params.acq.nacquisitions
                  for partition = 1:params.acq.npartitions
                    counter = counter+1;    
%                    waitbar(counter/total_readouts,h)
                    
                    params.phase=cardiac_phase;
                   % temp=read_n4_data_slice_ver5(rawfilename, params, slice, cardiac_phase, repetition, dataset, echo1, acq, 1);
                    temp=read_n4_data_slice_ver5(rawfilename, params, slice, cardiac_phase, repetition, dataset, echo1, acq, partition);% Alan adds partition
                    mrdata.raw(:,:,partition,:,cardiac_phase,dataset,acq,echo1,repetition,slice) = temp;
                    %tmp=fftshift(fft2(window.*cout));
                    %magdiff=rss(cout(:,:,:,cardiac_phase,1)-cout(:,:,:,cardiac_phase,2),3);
                    %p=squeeze(angle(sum(conj(cout(:,:,:,cardiac_phase,1)).*cout(:,:,:,cardiac_phase,2),3)));
                    %imagescn(cat(3,I,magdiff,p),[],[1 4],9) ; colormap(gray); drawnow
                  end
                end
              end
            end
          end 
        end % slice loop
        tic
    end
%    close(h)
else
	disp(['Can not read the file ' FILENAMEfull]);
	data = [];
end

mrdata.params = params;
mrdata.readoutindicator = readoutindicator;

return
      






