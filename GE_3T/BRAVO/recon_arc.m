function recon_arc(pfilePath, kacq_file, outputDir)

    if ~ exist('outputDir','var') || isempty(outputDir)
        outputDir = pwd;
    end


    % Load Pfile
    clear GERecon
    pfile = GERecon('Pfile.Load', pfilePath);
    GERecon('Pfile.SetActive',pfile);
    header = GERecon('Pfile.Header', pfile);


    % if length(varargin) == 1
    %     % Load Arc Sampling Pattern (kacq_yz.txt)
    %     GERecon('Arc.LoadKacq', varargin{1});
    % end

    % Load KSpace. Since 3D Arc Pfiles contain space for the zipped
    % slices (even though the data is irrelevant), only pull out
    % the true acquired K-Space. Z-transform will zip the slices
    % out to the expected extent.
    acquiredSlices = pfile.slicesPerPass / header.RawHeader.zip_factor;


    % 3D Scaling Factor
    scaleFactor = header.RawHeader.user0;
    if header.RawHeader.a3dscale > 0
        scaleFactor = scaleFactor * header.RawHeader.a3dscale;
    end

    scaleFactor = pfile.slicesPerPass / scaleFactor;



    % extract the kSpace
    % Synthesize KSpace to get full kSpace
    kSpace = zeros(pfile.xRes, pfile.yRes, acquiredSlices, pfile.channels, pfile.echoes, pfile.passes,'single');

    if exist('kacq_file','var') && ~isempty(kacq_file)
        GERecon('Arc.LoadKacq', kacq_file);
    end

    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
            for slice = 1:acquiredSlices            
                sliceInfo.pass = pass;
                sliceInfo.sliceInPass = slice;            
                for channel = 1:pfile.channels
                    % Load K-Space
                    kSpace(:,:,slice,channel,echo,pass) = GERecon('Pfile.KSpace', sliceInfo, echo, channel, pfile);
                end
            end
            kSpace(:,:,:,:,echo,pass) = GERecon('Arc.Synthesize', kSpace(:,:,:,:,echo,pass));
        end
    end



    % ASSET recon
    % change the p-file header of ASSET
    setenv('pfilePath',pfilePath);
    % unix('/Users/hongfusun/bin/orchestra/BuildOutputs/bin/HS_ModHeader --pfile $pfilePath');
    % unix('/Users/uqhsun8/bin/orchestra/BuildOutputs/bin/HS_ModHeader --pfile $pfilePath');
    unix('~/bin/orchestra/BuildOutputs/bin/HS_ModHeader --pfile $pfilePath');
    pfilePath=[pfilePath '.mod'];
    % Load Pfile
    clear GERecon
    pfile = GERecon('Pfile.Load', pfilePath);
    GERecon('Pfile.SetActive',pfile);
    header = GERecon('Pfile.Header', pfile);


    % image recon
    % Scale
    kSpace = kSpace * scaleFactor;
    % channelImages = zeros(pfile.xRes, pfile.yRes, acquiredSlices, pfile.channels, pfile.echoes, pfile.passes,'single');
    %%%%% 1 %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Default GE method %%%%%%%%%%%%%%%%%%%%%%%
    % Transform Across Slices
    kSpace = ifft(kSpace, pfile.slicesPerPass, 3);

    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
            for slice = 1:pfile.slicesPerPass
                for channel = 1:pfile.channels
                    % Transform K-Space
                    channelImages(:,:,channel) = GERecon('Transform', kSpace(:,:,slice,channel,echo,pass));
                end
                % Get slice information (corners and orientation) for this slice location
                sliceInfo.pass = pass;
                sliceInfo.sliceInPass = slice;
                %info = GERecon('Pfile.Info', sliceInfo);
                info = GERecon('Pfile.Info', slice);
%                 info = GERecon('Pfile.Corners', slice);
                % unaliasedImage(:,:,slice,echo,pass) = GERecon('Asset.Unalias', channelImages, info);
                unaliasedImage(:,:,:,slice,echo,pass) = channelImages;
            end
        end
    end


    clear kSpace
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % correct for phase chopping
    unaliasedImage = fft(fft(fft(fftshift(fftshift(fftshift(ifft(ifft(ifft(unaliasedImage,[],1),[],2),[],3),1),2),3),[],1),[],2),[],3);

    mkdir(outputDir);
    nii=make_nii(abs(unaliasedImage));
    save_nii(nii,[outputDir '/unaliasedImage_mag.nii']);
    nii=make_nii(angle(unaliasedImage));
    save_nii(nii,[outputDir '/unaliasedImage_pha.nii']);


    % save the matlab matrix for later use
    save([outputDir '/unaliasedImage'],'unaliasedImage','-v7.3');
