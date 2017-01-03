function [img_all, params] = epi15_recon(rawfile,ph_corr)



% read in raw kspace data
% [data, phascor1d, phascor2d, noise, patrefscan, ...
%     patrefscan_phascor, phasestabscan, refphasestabscan] ...
%         = read_meas_dat(filename);
[data, phascor1d] = read_meas_dat(rawfile);
datasize = size(data);
phascor1dsize = size(phascor1d);

data_all = data;
phascor1d_all = phascor1d;


% YAPS = read_meas_prot(filename); % NOT working!
% still use siem_read to fetch parameters
[pathstr, name, ext] = fileparts(rawfile);
mrdata = siem_read_params({[pathstr,filesep],[name,ext]});
params = mrdata.params;
sKSpace = params.protocol_header.sKSpace;



for timeSeries =1:datasize(7)
    data = permute(reshape(data_all(:,:,:,:,:,:,timeSeries,:,:,:),[datasize(1), datasize(2), datasize(3), datasize(8), datasize(10)]),[2 1 5 3 4]);
    k = sum(data,5); % PC x RO x NS x RV

    phascor1d = permute(reshape(phascor1d_all(:,:,:,:,:,:,timeSeries,:,:,:),[phascor1dsize(1), phascor1dsize(2), phascor1dsize(3), phascor1dsize(8), phascor1dsize(10)]),[2 1 5 3 4]);
    ref = sum(phascor1d,5); % 3 PC lines of ref scan




% phase correction (N/2 ghost)
switch ph_corr

    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % linear phase correction
        % find the peak of each readout line
        [PC, RO, NS, RV] = size(k);

        ref_f = fftshift(fft(fftshift(ref,2),[],2),2);
        odd_ref_f = (ref_f(1,:,:,:) + ref_f(3,:,:,:))/2;
        odd_ref = squeeze(ifftshift(ifft(ifftshift(odd_ref_f,2),[],2),2));

        [C,J] =  max(abs(odd_ref));
        I = zeros(1,NS,RV);
        A = zeros(1,NS,RV);

        % parabola fit 
        for i = 1:NS
            for j = 1:RV
                P = polyfit(J(1,i,j)-2:J(1,i,j)+2,abs(odd_ref(J(1,i,j)-2:J(1,i,j)+2,i,j))',2);
                I(1,i,j) = -P(2)/(2*P(1));
                % Q = polyfit(J(1,i,j)-2:J(1,i,j)+2,angle(odd_ref(J(1,i,j)-2:J(1,i,j)+2,i,j))',1);
                % A(1,i,j) = Q(1)*I(1,i,j)+Q(2);
            end
        end

        IS = RO/2 + 1 - I; % vox shifted

        % F_odd = exp(-1i*repmat(A,[RO 1 1])).*...
            % exp(-1i*2*pi*repmat(IS,[RO 1 1]).*repmat((-1/2:1/RO:1/2-1/RO)',[1 NS RV]));
        F_odd = exp(-1i*2*pi*repmat(IS,[RO 1 1]).*repmat((-1/2:1/RO:1/2-1/RO)',[1 NS RV]));


        even_ref = squeeze(ref(2,:,:,:));
        [C,J] = max(abs(even_ref));

        for i = 1:NS
            for j = 1:RV
                P = polyfit(J(1,i,j)-2:J(1,i,j)+2,abs(even_ref(J(1,i,j)-2:J(1,i,j)+2,i,j))',2);
                I(1,i,j) = -P(2)/(2*P(1));
                % Q = polyfit(J(1,i,j)-2:J(1,i,j)+2,angle(even_ref(J(1,i,j)-2:J(1,i,j)+2,i,j))',1);
                % A(1,i,j) = Q(1)*I(1,i,j)+Q(2);
            end
        end


        IS = RO/2 + 1 - I; % vox shifted

        % F_even = exp(-1i*repmat(A,[RO 1 1])).*...
            % exp(-1i*2*pi*repmat(IS,[RO 1 1]).*repmat((-1/2:1/RO:1/2-1/RO)',[1 NS RV]));
        F_even = exp(-1i*2*pi*repmat(IS,[RO 1 1]).*repmat((-1/2:1/RO:1/2-1/RO)',[1 NS RV]));


        %   Apply phase shift to data in hybrid space
        k = fftshift(fft(fftshift(k,2),[],2),2);
        for i = 2:2:datasize(2)
            k(i,:,:,:) = k(i,:,:,:) .* reshape(F_odd,[1 RO NS RV]);
        end
        for i = 1:2:datasize(2)
            k(i,:,:,:) = k(i,:,:,:) .* reshape(F_even,[1 RO NS RV]);
        end
        k = ifftshift(ifft(ifftshift(k,2),[],2),2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % non-linear phase correction
        k = fftshift(fft(fftshift(k,2),[],2),2);

        ref_f = fftshift(fft(fftshift(ref,2),[],2),2);
        odd_ref_f = ref_f(1,:,:,:) + ref_f(3,:,:,:);
        even_ref_f = ref_f(2,:,:,:);

        for i = 2:2:datasize(2) % start from the bottom
            k(i,:,:,:) = k(i,:,:,:).*exp(-1i*angle(odd_ref_f));
        end
        for i = 1:2:datasize(2)
            k(i,:,:,:) = k(i,:,:,:).*exp(-1i*angle(even_ref_f));
        end
        k = ifftshift(ifft(ifftshift(k,2),[],2),2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % linear MTF correction
        ref_f = fftshift(fft(fftshift(ref,2),[],2),2);
        temp1 = ref_f(1,:,:,:)./ref_f(2,:,:,:);
        temp2 = ref_f(3,:,:,:)./ref_f(2,:,:,:);
        phrmp = angle(temp1+temp2);

        % linear fit of phrmp (only ROI)
        % find the edge of signal
        std_ref = zeros(1,datasize(1)-10,datasize(10),datasize(3));
        ROI_a = zeros(datasize(10),datasize(3));
        ROI_b = ROI_a;

        diff = phrmp(:,[2:datasize(1),end],:,:) - phrmp;

        for m = 1:datasize(3)
            for j = 1:datasize(10)
                for i = 1:datasize(1)-10
                    std_ref(1,i,j,m) = std(diff(1,i:i+10,j,m));
                end
                
                ROI_a(j,m) = find(std_ref(1,:,j,m)<0.1,1,'first');
                ROI_b(j,m) = find(std_ref(1,:,j,m)<0.1,1,'last') + 10;
            end
        end

        phrmp_fit = zeros(1,datasize(1),datasize(10),datasize(3));
        % linear fit phrmp
        for j = 1:datasize(10)
            for m = 1:datasize(3)
                p = polyfit(ROI_a(j,m):ROI_b(j,m), phrmp(1,ROI_a(j,m):ROI_b(j,m),j,m), 1);
                phrmp_fit(1,:,j,m) = p(1)*(1:datasize(1)) + p(2);
            end
        end


        kr = fftshift(fft(fftshift(k,2),[],2),2);
        for i = 1:2:datasize(2)
            kr(i,:,:,:) = kr(i,:,:,:).*exp(1i*phrmp_fit/2);
            % kr(i,:,:,:) = kr(i,:,:,:).*exp(1i*phrmp/2);

        end
        for i = 2:2:datasize(2)
            kr(i,:,:,:) = kr(i,:,:,:).*exp(-1i*phrmp_fit/2);
            % kr(i,:,:,:) = kr(i,:,:,:).*exp(-1i*phrmp/2);
        end

        kr = ifftshift(ifft(ifftshift(kr,2),[],2),2);
        k = kr;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



% regridding the k-space in readout
if (params.protocol_header.m_alRegridMode(1) == 2) % 2 = REGRID_TRAPEZOIDAL

    [PC, RO, NS, RV] = size(k);
    k = reshape(permute(k,[2 1 3 4]),RO,[]); % RO x []

    RampupTime = params.protocol_header.m_alRegridRampupTime(1);
    FlattopTime = params.protocol_header.m_alRegridFlattopTime(1);
    RampdownTime = params.protocol_header.m_alRegridRampdownTime(1);
    DelaySamplesTime = params.protocol_header.m_alRegridDelaySamplesTime(1);
    ADCDuration = params.protocol_header.m_aflRegridADCDuration(1);
    DestSamples = params.protocol_header.m_alRegridDestSamples(1);


    % build the trapezoidal waveform
    rotrap = ones(1,RampupTime+FlattopTime+RampdownTime,'single');
    roramp = single(0:1/(RampupTime-1):1);
    rotrap(1:RampupTime)= roramp;
    rotrap(RampupTime+FlattopTime+1:end) = fliplr(roramp);

    % cut off the unused parts
    rotrap = rotrap(DelaySamplesTime+1:end);
    rotrap = rotrap(1:floor(ADCDuration)); % eja: floor added for VD11

    % integrate
    trapint = zeros(size(rotrap,2),1,'single');
    for z=1:size(rotrap,2)
        trapint(z) = sum(rotrap(1:z));
    end

    % assemble the desired k-space trajectory
    % add a point on the beginning and end since otherwise the
    %     interp1 function goes wacko
    destTraj = single(0:sum(rotrap)/(DestSamples+1):sum(rotrap));


    % interpolate 
    nSamples = size(k,1);
    actualDwell = ADCDuration/nSamples;
    actTraj = interp1(1:size(trapint,1),trapint,1:actualDwell:ADCDuration,'linear').';

    for i = 1:size(k,2)
        ctrc = k(:,i);

        filterf = (cumsum(destTraj(2:end-1)) - cumsum(actTraj'));
        filterf = filterf / sqrt(1/size(ctrc,1) * sum( abs(filterf).^2) );

        ctrc = interp1(actTraj,ctrc.*filterf',destTraj,'linear');
        % ctrc = interp1(actTraj,ctrc,destTraj.','linear');

        ctrc = ctrc(2:end-1);
        ctrc(isnan(ctrc)) = 0;
        k(:,i) = ctrc;

    end

    k = permute(reshape(k,[RO, PC, NS, RV]),[2 1 3 4]);
end



% partial fourier
sz = size(k);
pf = nan;
flag = sKSpace.ucPhasePartialFourier;
flag(isspace(flag))=[];
switch (flag)
    case {'0x1', 1}  % PF_HALF
        disp('Half fourier is unhandled :(')
    case {'0x2', 2}  % PF_5_8
        pf = 5/8;
    case {'0x4', 4}  % PF_6_8
        pf = 6/8;
    case {'0x8', 8}  % PF_7_8
        pf = 7/8;
    case {'0x10', 10} % PF_OFF
    case {'0x20', 20} % PF_AUTO
    case {'0x16', 16} % "none" in VD??
end
if ~isnan(pf)
    disp(sprintf('Partial Fourier: %d/8', pf*8));
    k = padarray(k, round(sz(1)*(1/(pf)-1)), 'pre');
end


% POCS
if size(k,4) > 1 && (~isnan(pf)) 
    [~, kspFull] = pocs(permute(k,[4 1 2 3]),20);
    k = permute(kspFull,[2 3 4 1]);
end


% phase resolution
if (sKSpace.dPhaseResolution ~= 1)
    disp(sprintf('Phase Resolution: %1.2f', sKSpace.dPhaseResolution));
    sz = size(k);
    pad_size = round(sz(1)*(1/sKSpace.dPhaseResolution-1));
    k = padarray(padarray(k,round(pad_size/2),'post'),pad_size-round(pad_size/2),'pre');
    % k = padarray(k, round((sz(1)-1)*(1/sKSpace.dPhaseResolution-1)/2));
end


% Asymmetric echo 
sz = size(k);
pad = sKSpace.lBaseResolution - sz(2)/2;
if pad
    disp(sprintf('Asymmetric echo: adding %d lines', pad*2));
    k = padarray(k, [0 pad*2], 'pre');
end


% 2D low-pass hann filter
% generate a 2d hamming low-pass filter
Nro = size(k,2);
Npe = size(k,1);
fw = 0.125;

x = hann(round(fw*Nro/2)*2);
x1 = [x(1:length(x)/2); ones([Nro-length(x),1]); x(length(x)/2+1:end)];
y = hann(round(fw*Npe/2)*2);
y1 = [y(1:length(y)/2); ones([Npe-length(y),1]); y(length(y)/2+1:end)];
        
[X,Y] = meshgrid(x1,y1);
Z = X.*Y;
Z = repmat(Z,[1 1 datasize(10) datasize(3)]);

k = k.*Z;



% fft to image space and remove oversampling
img = zeros(size(k));
for i = 1:size(k,3)*size(k,4)
    tmp = k(:,:,i);
    img(:,:,i) = fftshift(fft2(fftshift(tmp)));
end

img = img(:,size(k,2)/4+1:size(k,2)/4*3,:,:);


% permute the matrix first two dimensions (EPI acquisition)
img = permute(img,[2 1 3 4]);
% flip the readout and slice dimension to match scanner frame
img = flipdim(flipdim(img,1),3);





img_all(:,:,:,:,timeSeries) = img;
end