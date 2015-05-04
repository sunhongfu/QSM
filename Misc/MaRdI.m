classdef MaRdI
%MaRdI Ma(t)-R-dI(com)
%   Dicom into Matlab for Siemens MRI data
%
% Img = MaRdI( dataLoadDirectory )
%
%   Img contains fields
%
%       .img    (array of images - 3D if there are multiple DICOMs in directory)
%
%       .Hdr    (header of the first DICOM file read by dir( dataLoadDirectory ) ) 
% .......
%
% 2015
% topfer@ualberta.ca 
% =========================================================================
%   TODO
%
%   Saving FIELD as DICOM
%
%   function to ROTATE img based on Img.Hdr.ImageOrientationPatient

properties
    img ;
    Hdr ;
end

methods

function Img = MaRdI( dataLoadDirectory )

    if ( nargin == 0 )
        help(mfilename); 
        return; 
    end;
    
    listOfDicoms           = dir( [ dataLoadDirectory '/*.dcm'] );
    Img.Hdr                = dicominfo( [ dataLoadDirectory '/' listOfDicoms(1).name] ) ;
    
    Img.Hdr.NumberOfSlices = length(listOfDicoms) ;
    
    Img.img                = dicomread( [ dataLoadDirectory '/' listOfDicoms(1).name] )  ;

    for sliceIndex = 2 : Img.Hdr.NumberOfSlices
        Img.img(:,:,sliceIndex) = dicomread([ dataLoadDirectory '/' listOfDicoms(sliceIndex).name]) ;
    end
    
    Img = Img.scaleimgtophysical( ) ;   

end

function Img = scaleimgtophysical( Img, undoFlag )

    if (nargin < 2) 
        undoFlag = 0 ;
    end 
    
    if  ~myisfield( Img.Hdr, 'RescaleIntercept' ) 
        Img.Hdr.RescaleIntercept = 0 ;
    end
    
    if ~myisfield( Img.Hdr, 'RescaleSlope' )
        Img.Hdr.RescaleSlope = 1 ;
    end


    if undoFlag ~= -1
        %     assert( ~isempty(strfind( Img.Hdr.ImageType, '\P\' )) ... 
        %         ||  ~isempty(strfind( Img.Hdr.ImageType, 'Phase' )) % is pre-scaled phase (i.e. the img has been loaded and manipulated already (scaled by pi))
        %         ||  ~isempty(strfind( Img.Hdr.ImageType, 'Frequency' )) % is frequency
            
            Img.img = (Img.Hdr.RescaleSlope .* double( Img.img ) + Img.Hdr.RescaleIntercept)/double(2^Img.Hdr.BitsStored) ;

            if ~isempty( strfind( Img.Hdr.ImageType, '\P\' ) ) % is raw SIEMENS phase
                Img.img                             = pi*Img.img ; % scale to rad
                Img.Hdr.ImageType                   = 'PHASE\RAW\' ; 
            end
            
            Img.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none

            if ~isempty( strfind( Img.Hdr.ImageType, 'FREQUENCY\' ) )
                Img.Hdr.PixelComponentPhysicalUnits = '0005H' ; % i.e. Hz
            end

    else
            Img.Hdr.RescaleIntercept = min( Img.img(:) ) ;
        
            Img.img                  = Img.img - Img.Hdr.RescaleIntercept ;
        
            Img.Hdr.RescaleSlope     = max( Img.img(:) )/double(2^Img.Hdr.BitsStored ) ;
        
            Img.img                  = uint16( Img.img / Img.Hdr.RescaleSlope ) ;

    end



end

function Img = zeropad( Img, padSize, padDirection )

    gridSizeOriginalImg = size(Img.img) ;

    if (padSize(1) ~= 0) || (padSize(2) ~= 0) || ~strcmp(padDirection, 'post') 
        error( 'Feature not yet compatible. Fix it!' ) 
    end

    Img.img = padarray( Img.img, padSize, 0, padDirection ) ; 
    % Update header
    Img.Hdr.Rows                 = size(Img.img, 1) ;
    Img.Hdr.Columns              = size(Img.img, 2) ;       
    Img.Hdr.NumberOfSlices       = size(Img.img, 3) ;

    % if ~strcmp( padDirection, 'post' )
    %         if padSize(1) > 0 
    %             Img.Hdr.ImagePositionPatient(1) = ...
    %                 Img.Hdr.ImagePositionPatient(1) - Img.Hdr.PixelSpacing(1)*padSize(1) ) ;
    %         end
    %         if padSize() > 0 
    %             Img.Hdr.ImagePositionPatient(1) = ...
    %                 Img.Hdr.ImagePositionPatient(1) - Img.Hdr.PixelSpacing(1)*padSize(1) ) ;
    %         end
end

function Img = scalephasetofrequency( Img, undoFlag )
    %SCALEPHASETOFREQUENCY
    %
    % Converts unwrapped phase [units:rad] to field [units: Hz]
    % 
    %   Field = scalephasetofrequency( UnwrappedPhase )
    %
    % .....
    %
    % UnwrappedPhase.Hdr.EchoTime              [units: ms]

    scalingFactor = 1/( 2*pi*Img.Hdr.EchoTime*(1E-3)  ) ;
    
    if (nargin < 2) || (undoFlag ~= -1)
        assert( strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0000H' ) )

        Img.img       = scalingFactor * Img.img ;
        Img.Hdr.PixelComponentPhysicalUnits = '0005H' ; % i.e. Hz
    
    elseif (undoFlag == -1)
        assert( strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0005H' ) )

        Img.img       = (scalingFactor^-1 )* Img.img ;
        Img.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none
    end

end

function GYRO = getgyromagneticratio( Img )
    %GETGYROMAGNETICRATIO
    %
    % Examines .Hdr of MaRdI-type Img for .ImagedNucleus and returns gyromagnetic ratio in units of rad/s/T 
    %
    % .....
    %
    % Gyro = getgyromagneticratio( Img )

    switch Img.Hdr.ImagedNucleus 
        case '1H' 
            GYRO = 267.513E6 ; 
        otherwise
            error('Unknown nucleus. Do something useful: add case and corresponding gyromagnetic ratio.') ;
        end
end

function Phase = unwrapphase( Phase, Mag, PreludeOptions )
    %UNWRAPPHASE
    % Interfaces with FSL prelude
    
    assert( strcmp( Phase.Hdr.PixelComponentPhysicalUnits, '0000H' ), 'SCALEPHASE2RAD() beforehand.' ) ;
    
    if myisfield( Phase.Hdr, 'MaskingImage') && ~isempty(Phase.Hdr.MaskingImage)
        PreludeOptions.mask      = Phase.Hdr.MaskingImage ;
    end
    
    PreludeOptions.voxelSize = Phase.getvoxelsize() ;   

    Phase.img = prelude( Phase.img, Mag.img, PreludeOptions ) ;

    % update header
    Img.Hdr.ImageType        = 'PHASE\UNWRAPPED\' ; 
end

function voxelSize = getvoxelsize( Img )
    voxelSize = [ Img.Hdr.PixelSpacing(1) Img.Hdr.PixelSpacing(2) Img.Hdr.SpacingBetweenSlices ] ;
end

function fieldOfView = getfieldofview( Img )
    fieldOfView = [ Img.Hdr.PixelSpacing(1) * double( Img.Hdr.Rows - 1 ), ...
                    Img.Hdr.PixelSpacing(2) * double( Img.Hdr.Columns - 1 ), ...
                    Img.Hdr.SpacingBetweenSlices * double( Img.Hdr.NumberOfSlices - 1) ] ;
end

function gridSize = getgridsize( Img )
    gridSize = [ Img.Hdr.Rows, Img.Hdr.Columns, Img.Hdr.NumberOfSlices ] ;
end

function numberOfVoxels = getnumberofvoxels( Img )
    numberOfVoxels = prod( Img.getgridsize( ) ) ;
end

function [rowDirection, colDirection, sliceDirection] = getimgorientation( Img )
    
    rowDirection   = 0 ;
    colDirection   = 0 ;
    sliceDirection = 0 ;
    
    directions = [1 2 3] ;

    sliceDirection = find( (Img.Hdr.ImagePositionPatient == Img.Hdr.SliceLocation) ) ;

    for k = 1 : 3
        if (k ~= sliceDirection)
            if rowDirection == 0
                rowDirection = k ;
            else
                colDirection = k ;
            end
        end
    end

end

function [X,Y,Z] = getvoxelpositions( Img )
    %RCS = (DICOM) Reference Coordinate System
    fieldOfView = getfieldofview( Img ) ;

    [RCS_Z, RCS_X, RCS_Y] = ndgrid( [ Img.Hdr.ImagePositionPatient(3) : ...
                                      -Img.Hdr.PixelSpacing(1) : ...
                                      Img.Hdr.ImagePositionPatient(3) - fieldOfView(1) ], ...
                                        ...
                                    [ Img.Hdr.ImagePositionPatient(1) : ...
                                      Img.Hdr.PixelSpacing(2) : ...
                                      Img.Hdr.ImagePositionPatient(1) + fieldOfView(2) ], ...
                                        ...
                                    [ Img.Hdr.ImagePositionPatient(2) : ...
                                      Img.Hdr.SpacingBetweenSlices : ...
                                      Img.Hdr.ImagePositionPatient(2) + fieldOfView(3) ] )  ;
    X = RCS_Z ;
    Y = RCS_X ;
    Z = RCS_Y ;
end


function Img = resliceimg( Img, X_1, Y_1, Z_1, interpolationMethod ) 

    %%------ 
    % Reslice to new resolution
    if nargin < 5
        interpolationMethod = 'linear' ;
    end
    
    [X_0, Y_0, Z_0] = Img.getvoxelpositions( ) ;
    
    Img.img = interp3( Y_0, X_0, Z_0, Img.img, Y_1, X_1, Z_1, interpolationMethod ) ;

    % if new positions are outside the range of the original, interp3 replaces array entries with NaN
    Img.img( isnan( Img.img ) ) = 0 ; 

    % Update header
    Img.Hdr.PixelSpacing(1)      = abs(X_1(2,1,1) - X_1(1,1,1)) ;
    Img.Hdr.PixelSpacing(2)      = abs(Y_1(1,2,1) - Y_1(1,1,1)) ;
    Img.Hdr.SpacingBetweenSlices = abs(Z_1(1,1,2) - Z_1(1,1,1)) ;

    % KNOWLEDGE OF IMG ORIENTATION REQUIRED
    Img.Hdr.ImagePositionPatient( 3 ) = X_1(1) ; 
    Img.Hdr.ImagePositionPatient( 1 ) = Y_1(1) ;
    Img.Hdr.ImagePositionPatient( 2 ) = Z_1(1) ;

    Img.Hdr.Rows                 = size(Img.img, 1) ;
    Img.Hdr.Columns              = size(Img.img, 2) ;       
    Img.Hdr.NumberOfSlices       = size(Img.img, 3) ;

end


function [] = write( Img, dataSaveDirectory )
    %WRITE Ma(t)-R-dI(com)
    % 
    %.....
    %   Syntax
    %
    %   WRITE( dataSaveDirectory )
    %.....

    fprintf(['\n Writing img to: ' dataSaveDirectory ' ... \n'])

    [~,~,~] = mkdir( dataSaveDirectory ) ;

     Img.Hdr = rmfield(Img.Hdr, 'PixelComponentPhysicalUnits') ;

    if myisfield(Img.Hdr, 'MaskingImage')
        Img.Hdr = rmfield(Img.Hdr, 'MaskingImage') ;
    end

    Img = scaleimgtophysical( Img, -1 ) ;

    for sliceIndex = 1 : Img.Hdr.NumberOfSlices 

        sliceSuffix = '0000' ;
        sliceSuffix( end + 1 - length(num2str(sliceIndex) ) : end ) = num2str( sliceIndex ) ;
        sliceSuffix = ['-' sliceSuffix '.dcm'] ;

        %Img.Hdr.SliceLocation = ( sliceIndex - 1 ) * Img.Hdr.SpacingBetweenSlices ;

        dicomwrite( Img.img(:,:,sliceIndex) , ...
            [dataSaveDirectory '/' Img.Hdr.PatientName.FamilyName sliceSuffix], ...
            Img.Hdr, ...
            'WritePrivate', false, ...
            'UseMetadataBitDepths', true, ...
            'CreateMode', 'Copy' ) ; 
    end
end

end





end
