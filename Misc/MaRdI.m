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
    
    Img.Hdr.NumberOfSlices = uint16( length(listOfDicoms) ) ;
    
    Img.img                = dicomread( [ dataLoadDirectory '/' listOfDicoms(1).name] )  ;

    for sliceIndex = 2 : Img.Hdr.NumberOfSlices
        Img.img(:,:,sliceIndex) = dicomread([ dataLoadDirectory '/' listOfDicoms(sliceIndex).name]) ;
    end
    
    Img = Img.scaleimgtophysical( ) ;   

    if ~myisfield( Img.Hdr, 'SpacingBetweenSlices' ) 
        Img.Hdr.SpacingBetweenSlices = Img.Hdr.SliceThickness ;
    end

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
%ZEROPAD
% Img = ZEROPAD( Img, padSize, padDirection )
%
% padSize = [ nZerosRows nZerosColumns nZerosSlices ] 
%
% padDirection == 'post' || 'pre' || 'both'
%
% -------
%   See also PADARRAY()
    
    gridSizeOriginalImg = size(Img.img) ;

    Img.img = padarray( Img.img, padSize, 0, padDirection ) ; 
    
    % Update header
    Img.Hdr.Rows                 = size(Img.img, 1) ;
    Img.Hdr.Columns              = size(Img.img, 2) ;       
    Img.Hdr.NumberOfSlices       = size(Img.img, 3) ;

    if ~strcmp( padDirection, 'post' )
    % update image position 
    % (i.e. location in DICOM RCS of 1st voxel in data array (.img))
        
        % the following could be more efficient 
        % at the expensive of clarity
        voxelSize = Img.getvoxelsize() ;

        x0   = Img.Hdr.ImagePositionPatient(1) ;
        y0   = Img.Hdr.ImagePositionPatient(2) ;
        z0   = Img.Hdr.ImagePositionPatient(3) ;
        nr   = padSize(1)   ;
        nc   = padSize(2)   ;
        ns   = padSize(3)   ;
        dr   = voxelSize(1) ;
        dc   = voxelSize(2) ;
        ds   = voxelSize(3) ;
        rHat = Img.Hdr.ImageOrientationPatient(1:3) ; % unit vector row dir
        cHat = Img.Hdr.ImageOrientationPatient(4:6) ; % unit vector column dir
        sHat = cross( rHat, cHat ) ; % unit vector slice direction

        % -1 because the zeros are padded before ('pre') 1st element (origin)        
        dx = -1*( rHat(1)*dr*nr + cHat(1)*dc*nc + sHat(1)*ds*ns ) ; 
        dy = -1*( rHat(2)*dr*nr + cHat(2)*dc*nc + sHat(2)*ds*ns ) ;
        dz = -1*( rHat(3)*dr*nr + cHat(3)*dc*nc + sHat(3)*ds*ns ) ;

        x1 = x0 + dx ; 
        y1 = y0 + dy ;
        z1 = z0 + dz ;

        Img.Hdr.ImagePositionPatient = [x1 y1 z1] ;
    end

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

function Phase = unwrapphase( Phase, varargin )
    %UNWRAPPHASE
    % Interface to FSL prelude or to Abdul-Rahman's path-based phase unwrapper
    %
    %    .......................
    %    To call FSL Prelude   
    % 
    %       Phase = UNWRAPPHASE( Phase, Mag, PreludeOptions )
    %       
    %       Phase and Mag are objects of type MaRdI
    %       See HELP PRELUDE for description of PreludeOptions 
    %   
    %    .......................
    %    To call the Abdul-Rahman unwrapper  
    % 
    %	    Phase = UNWRAPPHASE( Phase, Options )
    %       Phase is an object of type MaRdI
    %       See HELP UNWRAP3D for description of Options 
    %
    %    .......................
    % 
    % Abdul-Rahman H, et al. Fast and robust three-dimensional best path phase
    % unwrapping algorithm. Applied Optics, Vol. 46, No. 26, pp. 6623-6635, 2007.

    assert( strcmp( Phase.Hdr.PixelComponentPhysicalUnits, '0000H' ), 'SCALEPHASE2RAD() beforehand.' ) ;
    
    isCallingPrelude = false ;
    Options.dummy    = [];
     
    if nargin == 1
        varargin{2} = Options ;
    end

   if nargin < 3 
    
        assert( myisfield( Phase.Hdr, 'MaskingImage') && ~isempty(Phase.Hdr.MaskingImage), ...
            'Logical masking array must be defined in Phase.Hdr.MaskingImage ' ) ;

        if nargin == 2
            Options = varargin{2} ;
        end
        Phase.img = unwrap3d( Phase.img, logical(Phase.Hdr.MaskingImage), Options ) ;

    elseif nargin == 3 
        isCallingPrelude = true ;

        Mag     = varargin{1} ;
        Options = varargin{2} ;

        if myisfield( Phase.Hdr, 'MaskingImage') && ~isempty(Phase.Hdr.MaskingImage)
            Options.mask = single( Phase.Hdr.MaskingImage ) ;
        end
        
        Options.voxelSize = Phase.getvoxelsize() ;   

        Phase.img = prelude( Phase.img, Mag.img, Options ) ;

    end

    % update header
    Img.Hdr.ImageType        = 'PHASE\UNWRAPPED\' ; 
end

function voxelSize = getvoxelsize( Img )
    voxelSize = [ Img.Hdr.PixelSpacing(1) Img.Hdr.PixelSpacing(2) Img.Hdr.SpacingBetweenSlices ] ;
end

% function Img = getfielddifferencemap( Img1, Img2 )
% % ???
%
%         assert( strcmp(Img1.Hdr.PixelComponentPhysicalUnits, '0005H') &&...
%                 strcmp(Img2.Hdr.PixelComponentPhysicalUnits, '0005H') ) ; % i.e. Hz
%
%         Img = Img2 ;
%         Img.Hdr.EchoTime = abs( Img2.Hdr.EchoTime - Img1.Hdr.EchoTime ) ;
%         
%         if Img2.Hdr.EchoTime > Img1.Hdr.EchoTime
%             Img.img = Img2.img - Img1.img ;
%         elseif Img2.Hdr.EchoTime < Img1.Hdr.EchoTime
%             Img.img = Img1.img - Img2.img ;
%         else
%             error('Echo times equal. Field difference should be zero.');
%         end
% end

function fieldOfView = getfieldofview( Img )
    fieldOfView = [ Img.Hdr.PixelSpacing(1) * double( Img.Hdr.Rows -1 ), ...
                    Img.Hdr.PixelSpacing(2) * double( Img.Hdr.Columns -1 ), ...
                    Img.Hdr.SpacingBetweenSlices * double( Img.Hdr.NumberOfSlices -1 ) ] ;
end

function gridSize = getgridsize( Img )
    gridSize = double( [ Img.Hdr.Rows, Img.Hdr.Columns, Img.Hdr.NumberOfSlices ] ) ;
end

function nVoxels = getnumberofvoxels( Img )
    nVoxels = prod( Img.getgridsize( ) ) ;
end

function [X,Y,Z] = getvoxelpositions( Img )
% GETVOXELPOSITIONS
% 
% [X,Y,Z] = GETVOXELPOSITIONS( Img ) 
%   Returns 3 3d arrays of doubles w each element containing the
%   RCS location ('Reference Coordinate System', DICOM standard)
%   of the corresponding voxel [units: mm]. 
    
    fieldOfView = getfieldofview( Img ) ;

    % [RCS_Z, RCS_X, RCS_Y] = ndgrid( [ Img.Hdr.ImagePositionPatient(3) : ...
    %                                   -Img.Hdr.PixelSpacing(1) : ...
    %                                   Img.Hdr.ImagePositionPatient(3) - fieldOfView(1) ], ...
    %                                     ...
    %                                 [ Img.Hdr.ImagePositionPatient(1) : ...
    %                                   Img.Hdr.PixelSpacing(2) : ...
    %                                   Img.Hdr.ImagePositionPatient(1) + fieldOfView(2) ], ...
    %                                     ...
    %                                 [ Img.Hdr.ImagePositionPatient(2) : ...
    %                                   Img.Hdr.SpacingBetweenSlices : ...
    %                                   Img.Hdr.ImagePositionPatient(2) + fieldOfView(3) ] )  ;
    % X = RCS_Z ;
    % Y = RCS_X ;
    % Z = RCS_Y ;

    voxelSize = Img.getvoxelsize() ;
    
    rHat = Img.Hdr.ImageOrientationPatient(1:3) ; % unit vector row dir
    cHat = Img.Hdr.ImageOrientationPatient(4:6) ; % unit vector column dir
    sHat = cross( cHat, rHat ) ; % unit vector slice direction
    
    % form DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
    %
    % If Anatomical Orientation Type (0010,2210) is absent or has a value of
    % BIPED, the x-axis is increasing to the left hand side of the patient. The
    % y-axis is increasing to the posterior side of the patient. The z-axis is
    % increasing toward the head of the patient.
    %
    % Arrays containing row, column, and slice indices of each voxel
    assert( ~myisfield( Img.Hdr, 'AnatomicalOrientationType' ) || ...
            strcmp( Img.Hdr.AnatomicalOrientationType, 'BIPED' ), ...
            'Error: AnatomicalOrientationType not supported.' ) ;
            
    [R,C,S] = ndgrid( [0:1:Img.Hdr.Rows-1], ...
                      [0:1:Img.Hdr.Columns-1], ...
                      [0:1:Img.Hdr.NumberOfSlices-1] ) ;

    %-------
    % SCALE to physical by sample spacing 
    % (i.e. grid size, i.e. effective voxel size) 
    R = voxelSize(1)*double(R);% voxelSize(1): spacing btw adjacent rows
    C = voxelSize(2)*double(C);% voxelSize(2): spacing btw adjacent columns
    S = voxelSize(3)*double(S);

    %-------
    % ROTATE to align row direction with x-axis, 
    % column direction with y-axis, slice with z-axis
    X1 = cHat(1)*R + rHat(1)*C + sHat(1)*S;
    Y1 = cHat(2)*R + rHat(2)*C + sHat(2)*S;
    Z1 = cHat(3)*R + rHat(3)*C + sHat(3)*S;

    %-------
    % TRANSLATE w.r.t. origin 
    % (i.e. location of 1st element: .img(1,1,1))
    X = Img.Hdr.ImagePositionPatient(1) + X1 ; 
    Y = Img.Hdr.ImagePositionPatient(2) + Y1 ; 
    Z = Img.Hdr.ImagePositionPatient(3) + Z1 ; 

end


function Img = resliceimg( Img, X_1, Y_1, Z_1, interpolationMethod ) 

    %%------ 
    % Reslice to new resolution
    if nargin < 5
        interpolationMethod = 'linear' ;
    end
    
    [X_0, Y_0, Z_0] = Img.getvoxelpositions( ) ;
    
    % Img.img = interp3( Y_0, X_0, Z_0, Img.img, Y_1, X_1, Z_1, interpolationMethod ) ;
    % Img.img = griddata( Y_0, X_0, Z_0, Img.img, Y_1, X_1, Z_1, interpolationMethod ) ;

    Img.img = griddata( Y_0, Z_0, X_0, Img.img, Y_1, Z_1, X_1, interpolationMethod ) ;
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
