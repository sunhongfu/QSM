function recon_arc_asset_slow(pfilePath, calibrationPfile, outputDir)

    if ~ exist('outputDir','var') || isempty(outputDir)
        outputDir = pwd;
    end
    % % mac
    % cd '/Users/hongfusun/DATA/p-files/oct7/raw/'
    % pfilePath='/Users/hongfusun/DATA/p-files/oct7/raw/P31232.7';
    % calibrationPfile='/Users/hongfusun/DATA/p-files/oct7/raw/P32328.7';

    % pfilePath='/Users/hongfusun/DATA/p-files/oct12/P44544.7';
    % pfilePath='/Users/hongfusun/DATA/p-files/oct12/P32256.7';

    % linux
    % pfilePath='/media/data/p-files/oct12/t1/P44544.7'

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
        end
    end


    % Synthesize KSpace to get full kSpace
    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
            kSpace(:,:,:,:,echo,pass) = GERecon('Arc.Synthesize', kSpace(:,:,:,:,echo,pass));       
        end
    end


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
                    channelImages(:,:,slice,channel,echo,pass) = GERecon('Transform', kSpace(:,:,slice,channel,echo,pass));
                end
            end
        end
    end

    nii=make_nii(abs(channelImages));
    save_nii(nii,'channelImage_mag.nii');

    nii=make_nii(angle(channelImages));
    save_nii(nii,'channelImage_pha.nii');

    clear kSpace
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % %%%%%% 2 %%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% IFFT method %%%%%%%%%%%%%%%%%%%%%%%
    % kSpace = padarray(kSpace,[18, 18], 'both');

    % channelImages = ifft(ifft(ifft(kSpace,[],1),[],2),[],3);
    % channelImages = fftshift(fftshift(channelImages,1),2);
    % channelImages = channelImages*256*sqrt(40);

    % nii=make_nii(abs(channelImages(:,:,:,:,4)));
    % save_nii(nii,'channelImage_mag_e4.nii');

    % nii=make_nii(angle(channelImages(:,:,:,:,4)));
    % save_nii(nii,'channelImage_pha_e4.nii');

    % clear kSpace
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    % %%%%%% 3 %%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Hongfu method %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % kSpace = padarray(kSpace,[18, 18], 'both');
    % kSpace = fftshift(kSpace,1);
    % kSpace = fftshift(kSpace,2);
    % kSpace = fftshift(kSpace,3);

    % channelImages = fft(fft(fft(kSpace,[],1),[],2),[],3);
    % channelImages = fftshift(fftshift(channelImages,1),2);
    % channelImages = channelImages/256/sqrt(40);

    % %%% seems like need flip to be the same as method 1 and 3, also circshift
    % channelImages = flipdim(flipdim(flipdim(channelImages,1),2),3);
    % channelImages = circshift(channelImages,[1 1 1]);

    % nii=make_nii(abs(channelImages(:,:,:,:,4)));
    % save_nii(nii,'channelImage_mag_e4.nii');

    % nii=make_nii(angle(channelImages(:,:,:,:,4)));
    % save_nii(nii,'channelImage_pha_e4.nii');

    % clear kSpace
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % %%%%%% 4 %%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Hongfu method no phase checkerboard %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % kSpace = fftshift(kSpace,1);
    % kSpace = fftshift(kSpace,2);
    % kSpace = fftshift(kSpace,3);

    % channelImages = fft(fft(fft(kSpace,[],1),[],2),[],3);
    % channelImages = fftshift(fftshift(channelImages,1),2);
    % channelImages = channelImages/256/sqrt(40);

    % % %%% seems like need flip to be the same as method 1 and 3, also circshift
    % % channelImages = flipdim(flipdim(flipdim(channelImages,1),2),3);
    % % channelImages = circshift(channelImages,[1 1 1]);

    % nii=make_nii(single(abs(channelImages)));
    % save_nii(nii,'channelImage_mag.nii');

    % nii=make_nii(single(angle(channelImages)));
    % save_nii(nii,'channelImage_pha.nii');

    % clear kSpace
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % % SVD coil combination
    % img_cmb = zeros(pfile.xRes, pfile.yRes, acquiredSlices);
    % for slice = 1:pfile.slicesPerPass
    %     [img_cmb(:,:,slice)] = adaptive_cmb_2d(squeeze(channelImages(:,:,slice,:)));
    % end

    % nii=make_nii(abs(img_cmb));
    % save_nii(nii,'combinedImage_mag.nii');

    % nii=make_nii(angle(img_cmb));
    % save_nii(nii,'combinedImage_pha.nii');




    % % combine magnitude images with SoS
    % % Apply Channel Combination
    % imsize = size(channelImages);
    % combinedImage = zeros([imsize(1:3), pfile.echoes, pfile.passes]);
    % finalImage = combinedImage;
    % for pass = 1:pfile.passes
    %     for echo = 1:pfile.echoes
    %         for slice = 1:acquiredSlices    
    %             % Get slice information (corners and orientation) for this slice location
    %             sliceInfo.pass = pass;
    %             sliceInfo.sliceInPass = slice;
    %             info = GERecon('Pfile.Info', sliceInfo);

    %             combinedImage(:,:,slice,echo,pass) = GERecon('SumOfSquares', squeeze(channelImages(:,:,slice,:,echo,pass)));

    %             % Create Magnitude Image
    %             magnitudeImage = abs(combinedImage(:,:,slice,echo,pass));

    %             % Apply Gradwarp
    %             gradwarpedImage = GERecon('Gradwarp', magnitudeImage, info.Corners);

    %             % Orient the image
    %             finalImage(:,:,slice,echo,pass) = GERecon('Orient', gradwarpedImage, info.Orientation);
    %         end
    %     end
    % end

    % nii=make_nii(abs(combinedImage(:,:,:,4)));
    % save_nii(nii,'combinedImage_mag_e4.nii');

    % nii=make_nii(abs(finalImage(:,:,:,4)));
    % save_nii(nii,'finalImage_mag_e4.nii');




    % ASSET recon
    % change the p-file header of ASSET
    setenv('pfilePath',pfilePath);
    unix('/Users/uqhsun8/bin/HS_ModHeader --pfile $pfilePath');
    pfilePath=[pfilePath '.mod'];
    % Load Pfile
    clear GERecon
    pfile = GERecon('Pfile.Load', pfilePath);
    GERecon('Pfile.SetActive',pfile);
    header = GERecon('Pfile.Header', pfile);


    % calibrationPfile
    GERecon('Calibration.Process', calibrationPfile);

    imsize = size(channelImages);

    unaliasedImage = zeros([imsize(1:3), pfile.echoes, pfile.passes]);
    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
            for slice = 1:acquiredSlices
                % Get slice information (corners and orientation) for this slice location
                sliceInfo.pass = pass;
                sliceInfo.sliceInPass = slice;
                info = GERecon('Pfile.Info', slice);
                unaliasedImage(:,:,slice,echo,pass) = GERecon('Asset.Unalias', squeeze(channelImages(:,:,slice,:,echo,pass)), info);
            end
        end
    end



    % correct for phase chopping
    unaliasedImage = fft(fft(fft(fftshift(fftshift(fftshift(ifft(ifft(ifft(unaliasedImage,[],1),[],2),[],3),1),2),3),[],1),[],2),[],3);

    nii=make_nii(abs(unaliasedImage));
    save_nii(nii,'unaliasedImage_mag.nii');
    nii=make_nii(angle(unaliasedImage));
    save_nii(nii,'unaliasedImage_pha.nii');


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

                % Save DICOMs
                imageNumber = ImageNumber(pass, info.Number, echo, pfile);
                filename = [outputDir '/DICOMs_real/realImage' num2str(imageNumber,'%03d') '.dcm'];
                GERecon('Dicom.Write', filename, finalRealImage, imageNumber, info.Orientation, info.Corners);
                filename = [outputDir '/DICOMs_imag/imagImage' num2str(imageNumber,'%03d') '.dcm'];
                GERecon('Dicom.Write', filename, finalImagImage, imageNumber, info.Orientation, info.Corners);
            end
        end
    end

end
