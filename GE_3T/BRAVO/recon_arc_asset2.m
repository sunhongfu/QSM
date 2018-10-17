function recon_arc_asset2(pfilePath, calibrationPfile, outputDir)

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
    unix('/Users/hongfusun/bin/orchestra-sdk-1.7-1/build/BuildOutputs/bin/HS_ModHeader --pfile $pfilePath');
    pfilePath=[pfilePath '.mod'];
    % Load Pfile
    clear GERecon
    pfile = GERecon('Pfile.Load', pfilePath);
    GERecon('Pfile.SetActive',pfile);
    header = GERecon('Pfile.Header', pfile);


    % calibrationPfile
    GERecon('Calibration.Process', calibrationPfile);




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
                info = GERecon('Pfile.Info', slice);
                unaliasedImage(:,:,slice,echo,pass) = GERecon('Asset.Unalias', channelImages, info);
            end
        end
    end


    clear kSpace
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % correct for phase chopping
    unaliasedImage = fft(fft(fft(fftshift(fftshift(fftshift(ifft(ifft(ifft(unaliasedImage,[],1),[],2),[],3),1),2),3),[],1),[],2),[],3);

    nii=make_nii(abs(unaliasedImage));
    save_nii(nii,'unaliasedImage_mag.nii');
    nii=make_nii(angle(unaliasedImage));
    save_nii(nii,'unaliasedImage_pha.nii');


    % save the matlab matrix for later use
    save('unaliasedImage','unaliasedImage');


    mkdir(outputDir);
    mkdir([outputDir '/DICOMs_real']);
    mkdir([outputDir '/DICOMs_imag']);
    % save DICOMs for QSM inputs
    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
            for slice = 1:acquiredSlices    
                % Get slice information (corners and orientation) for this slice location
                sliceInfo.pass = pass;
                sliceInfo.sliceInPass = slice;
                info = GERecon('Pfile.Info', sliceInfo);

                realImage = real(unaliasedImage(:,:,slice,echo,pass));
                imagImage = imag(unaliasedImage(:,:,slice,echo,pass));

                % Apply Gradwarp
                gradwarpedRealImage = GERecon('Gradwarp', realImage, info.Corners);
                gradwarpedImagImage = GERecon('Gradwarp', imagImage, info.Corners);

                % Orient the image
                finalRealImage = GERecon('Orient', gradwarpedRealImage, info.Orientation);
                finalImagImage = GERecon('Orient', gradwarpedImagImage, info.Orientation);

                tenum.Group = hex2dec('0018');
                tenum.Element = hex2dec('0086');
                tenum.VRType = 'IS';
                teval.Group = hex2dec('0018');
                teval.Element = hex2dec('0081');
                teval.VRType = 'DS';

                tenum.Value = num2str( echo );
                teval.Value = num2str( header.RawHeader.echotimes(echo) );

                % Save DICOMs
                imageNumber = ImageNumber(pass, info.Number, echo, pfile);
                filename = [outputDir '/DICOMs_real/realImage' num2str(imageNumber,'%03d') '.dcm'];
                GERecon('Dicom.Write', filename, finalRealImage, imageNumber, info.Orientation, info.Corners, (1000), 'desp', tenum, teval);
                filename = [outputDir '/DICOMs_imag/imagImage' num2str(imageNumber,'%03d') '.dcm'];
                GERecon('Dicom.Write', filename, finalImagImage, imageNumber, info.Orientation, info.Corners, (1000), 'desp', tenum, teval);
            end
        end
    end

end



function number = ImageNumber(pass, slice, echo, pfile)
% Image numbering scheme (P = Phase; S = Slice; E = Echo):
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn

    % Need to map the legacy "pass" number to a phase number
    numPassesPerPhase = fix(pfile.passes / pfile.phases);
    phase = fix(pass / numPassesPerPhase);

    slicesPerPhase = pfile.slicesPerPass * numPassesPerPhase * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1) + 1;
end


