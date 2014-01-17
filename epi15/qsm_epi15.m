function qsm_epi15(path_in)
% 1.5T EPI recon
init_dir = pwd;

PHASE_CORRECTION = 3;


if ~ exist('path_in','var') || isempty(path_in)
    path_in = pwd;
end

listing = dir([path_in '/*.out']);
pathstr = cd(cd(path_in));

for i = 1:size(listing,1)

    filename = listing(i).name;

    % read in raw kspace data
    % [data, phascor1d, phascor2d, noise, patrefscan, ...
    %     patrefscan_phascor, phasestabscan, refphasestabscan] ...
    %         = read_meas_dat(filename);
    [data, phascor1d] = read_meas_dat([pathstr,filesep,filename]);

    data = permute(squeeze(data),[2 1 5 3 4]);
    k = sum(data,5); % PC x RO x NS x RV

    phascor1d = permute(squeeze(phascor1d),[2 1 5 3 4]);
    ref = sum(phascor1d,5); % 3 PC lines of ref scan


    % YAPS = read_meas_prot(filename); % WHY NOT working?
    % still use siem_read to fetch parameters
    rawfile = {[pathstr,filesep],filename};

    mrdata = siem_read_v2(rawfile);
    params = mrdata.params;
    % k = mrdata.raw;
    % k = squeeze(k);
    % k = permute(k,[1 2 4 3]);


    % define directories
    path_out = [filename,'_QSM'];
    mkdir(pathstr,path_out);

    cd([pathstr,filesep,path_out]);

    % phase correction (N/2 ghost)
    switch PHASE_CORRECTION

        case 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % linear phase correction
            % find the peak of each readout line
            [PC, RO, NS, RV] = size(k);

            ref_f = fftshift(fft(fftshift(ref,2),[],2),2);
            odd_ref_f = (ref_f(1,:,:,:) + ref_f(3,:,:,:))/2;
            odd_ref = squeeze(ifftshift(ifft(ifftshift(odd_ref_f,2),[],2),2));


            % odd_ref = squeeze(ref(1,:,:,:)+ref(3,:,:,:));
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
            for i = 2:2:112
                k(i,:,:,:) = k(i,:,:,:) .* reshape(F_odd,[1 RO NS RV]);
            end
            for i = 1:2:112
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

            for i = 2:2:112 % start from the bottom
                k(i,:,:,:) = k(i,:,:,:).*exp(-1i*angle(odd_ref_f));
            end

            for i = 1:2:112
                k(i,:,:,:) = k(i,:,:,:).*exp(-1i*angle(even_ref_f));
            end

            k = ifftshift(ifft(ifftshift(k,2),[],2),2);


            % % additional MTF correction method
            % % self phase correction for ref scans
            % odd_ref_c = abs(odd_ref_f);
            % even_ref_c = abs(even_ref_f);


            % F_MTF = even_ref_c./odd_ref_c;

            % k = fftshift(fft(fftshift(k,2),[],2),2);

            % for i = 2:2:112
            %     k(i,:,:,:) = k(i,:,:,:).*F_MTF;
            % end

            % k = ifftshift(ifft(ifftshift(k,2),[],2),2);
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
            std_ref = zeros(1,246,60,8);
            ROI_a = zeros(60,8);
            ROI_b = ROI_a;

            diff = phrmp(:,[2:256,end],:,:) - phrmp;

            for m = 1:8
                for j = 1:60
                    for i = 1:246
                        std_ref(1,i,j,m) = std(diff(1,i:i+10,j,m));
                    end
                    
                    ROI_a(j,m) = find(std_ref(1,:,j,m)<0.1,1,'first');
                    ROI_b(j,m) = find(std_ref(1,:,j,m)<0.1,1,'last') + 10;
                end
            end

            phrmp_fit = zeros(1,256,60,8);
            % linear fit phrmp
            for j = 1:60
                for m = 1:8
                    p = polyfit(ROI_a(j,m):ROI_b(j,m), phrmp(1,ROI_a(j,m):ROI_b(j,m),j,m), 1);
                    phrmp_fit(1,:,j,m) = p(1)*(1:256) + p(2);
                end
            end


            kr = fftshift(fft(fftshift(k,2),[],2),2);
            for i = 1:2:112
                kr(i,:,:,:) = kr(i,:,:,:).*exp(1i*phrmp_fit/2);
                % kr(i,:,:,:) = kr(i,:,:,:).*exp(1i*phrmp/2);

            end
            for i = 2:2:112
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
    flag = params.protocol_header.sKSpace.ucPhasePartialFourier;
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
    [im, kspFull] = pocs(permute(k,[4 1 2 3]),20);
    k = permute(kspFull,[2 3 4 1]);


    % phase resolution
    if (params.protocol_header.sKSpace.dPhaseResolution ~= 1)
        disp(sprintf('Phase Resolution: %1.2f', params.protocol_header.sKSpace.dPhaseResolution));
        sz = size(k);
        k = padarray(k, round((sz(1)-1)*(1/params.protocol_header.sKSpace.dPhaseResolution-1)/2));
    end


    % Asymmetric echo 
    sz = size(k);
    pad = params.protocol_header.sKSpace.lBaseResolution - sz(2)/2;
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
    Z = repmat(Z,[1 1 60 8 ]);

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
    % flip the readout to match qsm_swi15 results
    img = flipdim(img,1);
    nii = make_nii(abs(img));
    save_nii(nii,'mag_all.nii');
    nii = make_nii(angle(img));
    save_nii(nii,'ph_all.nii');


    % size and resolution
    [Nro,Npe,Ns,~] = size(img);
    FOV = params.protocol_header.sSliceArray.asSlice{1};
    vox = [FOV.dReadoutFOV/Nro, FOV.dPhaseFOV/Npe,  FOV.dThickness];


    % % combine coils
    cref = 8; % reference coil
    radi = 5; % kernel size

    img_cmb = sense_se(img,vox,cref,radi);
    nii = make_nii(abs(img_cmb),vox);
    save_nii(nii,'mag.nii');
    nii = make_nii(angle(img_cmb),vox);
    save_nii(nii,'ph.nii');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % combine coils
    % % 
    % img_cmb = zeros(Nro,Npe,Ns);
    % matlabpool open
    % parfor i = 1:Ns
    %     img_cmb(:,:,i) = coilCombinePar(img(:,:,i,:));
    % end
    % matlabpool close
    % nii = make_nii(abs(img_cmb),vox);
    % save_nii(nii,'mag.nii');
    % nii = make_nii(angle(img_cmb),vox);
    % save_nii(nii,'ph.nii');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ! bet mag.nii BET -f 0.4 -m -R;
    ! gunzip -f BET.nii.gz;
    ! gunzip -f BET_mask.nii.gz;
    nii = load_nii('BET_mask.nii');
    mask = double(nii.img);
    img_cmb = img_cmb.*mask;

    ! prelude -a mag.nii -p ph.nii -u unph.nii -m BET_mask.nii -n 8;
    ! gunzip -f unph.nii.gz;
    nii = load_nii('unph.nii');
    unph = double(nii.img);


    % background field removal
    ker_rad = 5; % convolution kernel radius size (mm)
    tik_reg = 1e-3; % tikhonov regularization

    % % (1) PDF
    % theta = -acos(params.protocol_header.sSliceArray.asSlice{1}.sNormal.dTra);
    % [lfs,mask_ero] = pdf(tfs,mask,vox,ker_rad,abs(img_cmb),theta);
    % nii = make_nii(lfs,vox);
    % save_nii(nii,'lfs_pdf.nii');

    % (2) RESHARP
    [lph,mask_ero] = resharp(unph,mask,vox,ker_rad,tik_reg);
    nii = make_nii(mask_ero,vox);
    save_nii(nii,'mask_ero.nii');
    nii = make_nii(lph,vox);
    save_nii(nii,'lph.nii');


    % normalize to ppm unit
    TE = params.protocol_header.alTE{1}/1e6;
    B_0 = params.protocol_header.m_flMagneticFieldStrength;
    gamma = 2.675222e8;
    lfs = lph/(gamma*TE*B_0)*1e6; % unit ppm
    nii = make_nii(lfs,vox);
    save_nii(nii,'lfs.nii');

    
    % susceptibility inversion
    % account for oblique slicing (head tilted)
    theta = -acos(params.protocol_header.sSliceArray.asSlice{1}.sNormal.dTra);

    tv_reg = 5e-4; % total variation regularization

    sus = tvdi(lfs,mask_ero,vox,tv_reg,abs(img_cmb),theta);
    % sus_final = sus.*mask_ero;
    % nii = make_nii(sus_final,vox);
    nii = make_nii(sus,vox);
    save_nii(nii,'sus.nii');


    % debugging purpose
    % save all the variables in all.mat
    save all.mat;

end % end of the very first for loop

cd(init_dir);








% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % MEDI inversion
% iMag = abs(img_cmb);
% Mask = mask_ero;
% matrix_size = size(Mask);
% voxel_size = vox;
% merit = 0;
% smv = 0;
% radius = 5;
% data_weighting = 1;
% gradient_weighting = 1;
% pad = 0;
% matrix_size0 = 0;
% Debug_Mode = 'NoDebug';
% cg_max_iter = 100;
% cg_tol = 0.01;
% max_iter = 10;
% tol_norm_ratio = 0.1;
% data_weighting_mode = data_weighting;
% gradient_weighting_mode = gradient_weighting;
% grad = @cgrad;
% div = @cdiv;
% B0_dir = [0,0,1]';
% tempn = ones(size(Mask)); %%%%%%%%%%%%%%% need update
% D=dipole_kernel(matrix_size, voxel_size, B0_dir);
% m = dataterm_mask(data_weighting_mode, tempn, Mask);
% RDF = lfs; %%%%%%%%%%%%%%%%%%%%%% need unit update
% b0 = m.*exp(1i*RDF);
% wG = gradient_mask(gradient_weighting_mode, iMag, Mask, grad, voxel_size);
% iter=0;
% x = zeros(matrix_size);
% res_norm_ratio = Inf;
% cost_data_history = zeros(1,max_iter);
% cost_reg_history = zeros(1,max_iter);
% e=0.000001; %a very small number to avoid /0
% badpoint = zeros(matrix_size);

% lambda = 1000;
% while (res_norm_ratio>tol_norm_ratio)&&(iter<max_iter)
% %while iter<max_iter
% tic
%     iter=iter+1;
%     Vr = 1./sqrt(abs(wG.*grad(real(x),voxel_size)).^2+e);
%     w = m.*exp(1i*ifftn(D.*fftn(x)));
%     reg = @(dx) div(wG.*(Vr.*(wG.*grad(real(dx),voxel_size))),voxel_size);
%     fidelity = @(dx)2*lambda*real(ifftn(D.*fftn(conj(w).*w.*real(ifftn(D.*fftn(dx))))));

%     A =  @(dx) reg(dx) + fidelity(dx);       
%     b = reg(x) + 2*lambda*real(ifftn(D.*fftn( conj(w).*conj(1i).*(w-b0))));



%     dx = real(cgsolve(A, -b, cg_tol, cg_max_iter, 0));
%     res_norm_ratio = norm(dx(:))/norm(x(:));
%     x = x + dx;

%     wres=m.*exp(1i*(real(ifftn(D.*fftn(x))))) - b0;

%     cost_data_history(iter) = norm(wres(:),2);
%     cost=abs(wG.*grad(x));
%     cost_reg_history(iter) = sum(cost(:));

    
%     if merit
%         wres = wres - mean(wres(Mask(:)==1));
%         a = wres(Mask(:)==1);
%         factor = std(abs(a))*6;
%         wres = abs(wres)/factor;
%         wres(wres<1) = 1;
%         badpoint(wres>1)=1;
%         N_std(Mask==1) = N_std(Mask==1).*wres(Mask==1).^2;
%         tempn = N_std;
%         if (smv)
%                 tempn = sqrt(SMV(tempn.^2, matrix_size, voxel_size, radius)+tempn.^2);
%         end
%         m = dataterm_mask(data_weighting_mode, tempn, Mask);
%         b0 = m.*exp(1i*RDF);
%     end
    
% toc
% end

% nii = make_nii(x,vox);
% save_nii(nii,'x.nii');


