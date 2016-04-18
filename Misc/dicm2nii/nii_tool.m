function varargout = nii_tool(cmd, varargin)
% Basic function to create, load and save NIfTI file.
% 
% rst = nii_tool(cmd, para);
% 
% To list all command, type
%  nii_tool ?
% 
% To get help information for each command, include '?' in cmd, for example:
%  nii_tool init?
%  nii_tool('init?')
% 
% Here is a list of all command:
% 
% nii_tool('default', 'version', 1, 'rgb_dim', 1);
% nii = nii_tool('init', img);
% nii_tool('save', nii, filename, force_3D);
% hdr = nii_tool('hdr', filename);
% img = nii_tool('img', filename_or_hdr);
% ext = nii_tool('ext', filename_or_hdr);
% nii = nii_tool('load', filename_or_hdr);
% nii = nii_tool('cat3D', filenames);
% nii_tool('RGBStyle', 'afni');
% 
% Detail for each command is described below.
% 
% oldVal = nii_tool('default', 'version', 1, 'rgb_dim', 1);
% oldVal = nii_tool('default', struct('version', 1, 'rgb_dim', 1));
% 
% - Set/query default NIfTI version and/or rgb_dim. To check the setting, run
% nii_tool('default') without other input. The input for 'default' command can
% be either a struct with fields of 'version' and/or 'rgb_dim', or
% parameter/value pairs. See nii_tool('RGBstyle') for meaning of rgb_dim.
% 
% Note that the setting will be saved for future use. If one wants to change the
% settting temporarily, it is better to return the oldVal, and to restore it
% after done:
% 
%  oldVal = nii_tool('default', 'version', 2); % set version 2 as default
%  % 'init' and 'save' NIfTI using above version
%  nii_tool('default', oldVal); % restore default setting
% 
% The default version setting affects 'init' command only. If you 'load' a NIfTI
% file, modify it, and then 'save' it, the version will be the same as the
% original file, unless it is changed explicitly (see help for 'save' command).
% All 'load' command ('load', 'hdr', 'ext', 'img') will read any version
% correctly, regardless of version setting.
% 
% 
% nii = nii_tool('init', img, RGB_dim);
% 
% - Initialize nii struct based on img, normally 3D or 4D array. Most fields in
% the returned nii.hdr contain default values, and need to be updated based on
% dicom or other information. Important ones include pixdim and s/qform_code and
% related parameters.
% 
% The NIfTI datatype will depend on data type of img. Most Matlab data types are
% supported, including 8/16/32/64 bit signed and unsigned integers, single and
% double floating numbers. Single/double complex and logical array are also
% supported.
% 
% nii_tool returns img with the same data type as it is stored, while numeric
% values in hdr are in double. This also applies to the struct returned by any
% nii_tool loading sub-functions.
% 
% The optional third input is needed only if img contains RGB/RGBA data. It
% specifies which dimension in img is for RGB or RGBA. In other words, if a
% non-empty third input is provided, img will be interpreted as RGB or RGBA
% data.
% 
% Another way to signify RGB/RGBA data is to permute color dim to 8th-dim of img
% (RGB_dim of 8 can be omitted then). Since NIfTI img can have up to 7 dim,
% nii_tool chooses to store RGB/RGBA in 8th dim. Although this looks lengthy
% (4th to 7th dim are often all ones), nii_tool can deal with up to 7 dim
% without causing any confusion. This is why the returned nii.img will always
% store RGB in 8th dim.
% 
% 
% hdr = nii_tool('hdr', filename);
% 
% - Return hdr struct of the provided NIfTI file.
% 
% 
% img = nii_tool('img', filename_or_hdr);
% 
% - Return image data in a NIfTI file. The second input can be NIfTI file name,
% or hdr struct returned by nii_tool('hdr', filename).
% 
% 
% ext = nii_tool('ext', filename_or_hdr);
% 
% - Return NIfTI extension in a NIfTI file. The second input can be NIfTI file
% name, or hdr struct returned by nii_tool('hdr', filename). The returned ext
% will have field 'edata_decoded' if 'ecode' is of known type, such as dicom
% (2), text (4 or 6) or Matlab (40).
% 
% Here is an example to add data in myFile.mat as extension to nii struct, which
% can be from 'init' or 'load':
% 
%  fid = fopen('myFile.mat'); % open the MAT file
%  myEdata = fread(fid, inf, '*uint8'); % load all bytes as byte column
%  fclose(fid);
%  len = int32(numel(myEdata)); % number of bytes in int32
%  myEdata = [typecast(len, 'uint8')'; myEdata]; % include len in myEdata
%  nii.ext.ecode = 40; % 40 for Matlab extension
%  nii.ext.edata = myEdata; % myEdata must be uint8 array
% 
% nii_tool will take care of rest when you 'save' nii to a file.
% 
% In case a NIfTI ext causes problem (for example, some FSL builds have problem
% in reading NIfTI img with ecode>30), one can remove the ext easily:
% 
%  nii = nii_tool('load', 'file_with_ext.nii'); % load the file with ext
%  nii.ext = []; % or nii = rmfield(nii, 'ext'); % remove ext
%  nii_tool('save', nii, 'file_without_ext.nii'); % save it
%
% 
% nii = nii_tool('load', filename_or_hdr);
% 
% - Load NIfTI file into nii struct. The returned struct includes NIfTI 'hdr'
% and 'img', as well as 'ext' if the file contains NIfTI extension.
% 
% 
% nii_tool('save', nii, filename, force_3D);
% 
% - Save struct nii into filename. The format of the file is determined by the
% file extension, such as .img, .nii, .img.gz, .nii.gz etc. If filename is not
% provided, nii.hdr.file_name must contain a file name. Note that 'save' command
% always overwrites file in case of name conflict.
% 
% If filename has no extension, '.nii' will be used as default.
% 
% If the 4th input, force_3D, is true (default false), the output file will be
% 3D only, which means multiple volume data will be split into multiple files.
% This is the format SPM likes. You can use this command to convert 4D into 3D
% by 'load' a 4D file, then 'save' it as 3D files. The 3D file names will have
% 5-digit like '_00001' appended to indicate volume index.
% 
% The NIfTI version can be set by nii_tool('default'). One can override the
% default version by specifying it in nii.hdr.version. To convert between
% versions, load a NIfTI file, specify new version, and save it. For example:
% 
%  nii = nii_tool('load', 'file_nifti1.nii'); % load version 1 file
%  nii.hdr.version = 2; % force to NIfTI-2
%  nii_tool('save', nii, 'file_nifti2.nii'); % save as version 2 file
% 
% Following example shows how to change data type of a nii file:
%  nii = nii_tool('load', 'file_int16.nii'); % load int16 type file
%  nii.img = single(nii.img); % change data type to single/float32
%  nii_tool('save', nii, 'file_float.nii'); % nii_tool will take care of hdr
% 
% 
% nii = nii_tool('cat3D', files);
% 
% - Concatenate SPM 3D files into a 4D dataset. The input 'files' can be cellstr
% with file names, or char with wildcards (* or ?). If it is cellstr, the volume
% order in the 4D data corresponds to those files. If wildcards are used, the
% volume order is based on alphabetical order of file names.
% 
% Note that the files to be concatenated must have the same datatype, dim, voxel
% size, scaling slope and intercept, transformation matrix, etc. This is true if
% files are for the same series. 
% 
% Following example shows how to convert a series of 3D files into a 4D file:
% 
%  nii = nii_tool('cat3D', './data/fSubj2-0003*.nii'); % load files for series 3 
%  nii_tool('save', nii, './data/fSubj2-0003_4D.nii'); % save as a 4D file
% 
% 
% oldStyle = nii_tool('RGBStyle', 'afni');
% 
% - Set/query the method to save/load RGB or RGBA NIfTI file. The default method
% can be set by nii_tool('default', [version rgb_dim]), where rgb_dim can be 1,
% 3 or 4, as explained below.
% 
% The default style is 'afni' style (or 1), which is defined by NIfTI standard,
% but is not well supported by fsl or mricron.
% 
% If the second input is set to 'mricron' (or 3), nii_tool will save/load file
% using the old RGB fashion (dim 3 for RGB). This works for mricron at least
% till version 4 August 2014.
% 
% If the second input is set to 'fsl' (or 4), nii_tool will save/load RGB or
% RGBA layer into 4th dimension (and the file is not encoded as RGB data, but as
% normal file). This violates the NIfTI rule, but it seems it is the only way,
% for now (till fsl version 5.0.8), to work for fslview.
% 
% If no new style (second input) is provided, it means to query the current
% style (one of 'afni', 'mricron' and 'fsl').
% 
% Following shows how to convert between mricron style and fsl style:
% 
%  nii_tool('RGBStyle', 'mricron'); % later load/save uses mricron style
%  nii = nii_tool('load', 'mricronStyle.nii'); % load dim3-RGB file
%  nii_tool('RGBStyle', 'fsl'); % switch to fsl style
%  nii_tool('save', nii, 'fslRGB.nii'); % fsl can read it as RGB
% 
% Note that, if one wants to convert fsl style (non-RGB file by NIfTI standard)
% to other styles, an extra step is needed to change the RGB dim from 4th to 8th
% dim before 'save':
% 
%  nii = nii_tool('load', 'fslStyleFile.nii'); % no need to set 'RGBStyle'
%  nii.img = permute(nii.img, [1:3 5:8 4]); % force it to be RGB data
%  nii_tool('RGBStyle', 'mricron'); % switch to RGB mricron style
%  nii_tool('save', nii, 'mricronRGB.nii'); % now mricron can read it as RGB
% 
% Also note that the setting by nii_tool('RGBStyle') is effective only for
% current Matlab session. If one clears all, or starts a new Matlab session, the
% default style by nii_tool('default') will take effect.
 
% More information for NIfTI format:
% Official NIfTI website: http://nifti.nimh.nih.gov/
% Another excellent site: http://brainder.org/2012/09/23/the-nifti-file-format/

% History (yymmdd)
% 150109 Write it based on Jimmy Shen's NIfTI tool (xiangrui.li@gmail.com)
% 150202 Include renamed pigz files for Windows to trick Matlab Central
% 150203 Fix closeFile and deleteTmpFile order
% 150205 Add hdr.machine: needed for .img fopen
% 150208 Add 4th input for 'save', allowing to save SPM 3D files
% 150210 Add 'cat3D' to load SPM 3D files
% 150226 Assign all 8 char for 'magic' (version 2 needs it)
% 150321 swapbytes(nByte) for ecode=40 with big endian
% 150401 Add 'default' to set/query version and rgb_dim default setting
% 150514 read_ext: decode txt edata by dicm2nii.m
% 150517 fhandle: provide a way to use gunzipOS etc from outside

persistent C para; % C columns: name, length, format, value, offset
if isempty(C), [C, para] = niiHeader; end

if ~ischar(cmd)
    error('Provide a string command as the first input for nii_tool');
end
if any(cmd=='?'), subFuncHelp(mfilename, cmd); return; end

if strcmpi(cmd, 'init')
    if nargin<2, error('nii_tool(''%s'') needs second input', cmd); end
    for i = 1:size(C,1), nii.hdr.(C{i,1}) = C{i,4}; end
    nii.img = varargin{1};
    if numel(size(nii.img))>8
        error('NIfTI img can have up to 7 dimension');
    end
    if nargin>2
        i = varargin{2};
        if i<0 || i>8 || mod(i,1)>0, error('Invalid RGB_dim number'); end
        nii.img = permute(nii.img, [1:i-1 i+1:8 i]); % RGB to dim8
    end
    varargout{1} = img2datatype(nii, para); % set datatype etc 
    
elseif strcmpi(cmd, 'save')
    if nargin<2, error('nii_tool(''%s'') needs second input', cmd); end
    nii = varargin{1};
    if ~isstruct(nii) || ~isfield(nii, 'hdr') || ~isfield(nii, 'img') 
        error(['nii_tool(''save'') needs a struct from nii_tool(''init'')' ...
            ' or nii_tool(''load'') as the second input']);
    end
    
    % Check file name to save
    if nargin>2
        fname = varargin{2};
        if length(fname)<5 || ~ischar(fname)
            error('Invalid name for NIfTI file: %s', fname);
        end
    elseif isfield(nii.hdr, 'file_name')
        fname = nii.hdr.file_name;
    else
        error('Provide a valid file name as the third input');
    end
    if ~ispc && strncmp(fname, '~/', 2) % matlab may err with this abbrevation
        fname = [getenv('HOME') fname(2:end)];
    end
    [pth, fname, fext] = fileparts(fname);
    do_gzip = strcmpi(fext, '.gz');
    if do_gzip
        [~, fname, fext] = fileparts(fname); % get .nii .img .hdr
    end
    if isempty(fext), fext = '.nii'; end % default single .nii file
    fname = fullfile(pth, fname); % without file ext
    isNii = strcmpi(fext, '.nii'); % will use .img/.hdr if not .nii
    
    % Deal with NIfTI version and sizeof_hdr
    niiVer = para.version;
    if isfield(nii.hdr, 'version'), niiVer = nii.hdr.version; end
    if niiVer == 1
        nii.hdr.sizeof_hdr = 348; % in case it was loaded from other version
    elseif niiVer == 2
        nii.hdr.sizeof_hdr = 540; % version 2
    else 
        error('Unsupported NIfTI version: %g', niiVer);
    end
    
    if niiVer ~= para.version
        C0 = niiHeader(niiVer);
    else
        C0 = C;
    end
    
    % Update datatype/bitpix/dim in case nii.img is changed
    [nii, fmt] = img2datatype(nii, para);

    % This 'if' block: lazy implementation SPM: split to 3D files
    if nargin>3 && ~isempty(varargin{3}) && varargin{3} && nii.hdr.dim(5)>1
        if do_gzip, fext = [fext '.gz']; end
        nii0 = nii;
        for i = 1:nii.hdr.dim(5)
            fname0 = sprintf('%s_%05g%s', fname, i, fext);
            nii0.img = nii.img(:,:,:,i,:,:,:,:); % one vol
            if i==1 && isfield(nii, 'ext'), nii0.ext = nii.ext;
            elseif i==2 && isfield(nii0, 'ext'), nii0 = rmfield(nii0, 'ext'); 
            end
            nii_tool('save', nii0, fname0);
        end
        return;
    end
        
    % re-arrange img for special datatype: RGB/RGBA/Complex.
    % mis-use unused data_type to store rgb_dim in nifti FILE
    if any(nii.hdr.datatype == [128 511 2304]) % RGB or RGBA
        nii.hdr.data_type = sprintf('rgb_dim=%g', para.rgb_dim);
        if para.rgb_dim == 1 % AFNI style
            nii.img = permute(nii.img, [8 1:7]);
        elseif para.rgb_dim == 3 % mricron
            nii.img = permute(nii.img, [1 2 8 3:7]);
        elseif para.rgb_dim == 4 % for fslview
            nii.img = permute(nii.img, [1:3 8 4:7]); % violate nii rule
            dim = size(nii.img);
            if numel(dim)>6 % dim7 is not 1
                i = find(dim(5:7)==1, 1, 'last') + 4;
                nii.img = permute(nii.img, [1:i-1 i+1:8 i]);
            end
            nii = img2datatype(nii, para); % changed to non-RGB datatype
        end
    elseif any(nii.hdr.datatype == [32 1792]) % complex single/double
        nii.img = [real(nii.img(:))'; imag(nii.img(:))'];
    end
    
    % Check nii extension: update esize to x16
    nExt = 0; esize = 0;
    nii.hdr.extension = [0 0 0 0]; % no nii ext
    if isfield(nii, 'ext') && isstruct(nii.ext) ...
            && isfield(nii.ext(1), 'edata') && ~isempty(nii.ext(1).edata)
        nExt = length(nii.ext);
        nii.hdr.extension = [1 0 0 0]; % there is nii ext
        for i = 1:nExt
            if ~isfield(nii.ext(i), 'ecode') || ~isfield(nii.ext(i), 'edata')
                error('NIfTI header ext struct must have ecode and edata');
            end
            
            n0 = numel(nii.ext(i).edata) + 8; % 8 byte for esize and ecode
            n1 = ceil(n0/16) * 16; % esize: multiple of 16
            nii.ext(i).esize = n1;
            nii.ext(i).edata(end+(1:n1-n0)) = 0; % pad zeros
            esize = esize + n1;
        end
    end
    
    % Set magic, vox_offset, and open file for .nii or .hdr
    if isNii
        % version 1 will take only the first 4
        nii.hdr.magic = sprintf('n+%g%s', niiVer, char([0 13 10 26 10]));
        nii.hdr.vox_offset = nii.hdr.sizeof_hdr + 4 + esize;
        fid = fopen([fname fext], 'w');
    else
        nii.hdr.magic = sprintf('ni%g%s', niiVer, char([0 13 10 26 10]));
        nii.hdr.vox_offset = 0;
        fid = fopen([fname '.hdr'], 'w');
    end
    
    % Write nii hdr
    for i = 1:size(C0,1)
        if isfield(nii.hdr, C0{i,1})
            val = nii.hdr.(C0{i,1});
        else % niiVer=2 omit some fields, also take care of other cases
            val = C0{i,4};
        end
        n = numel(val);
        len = C0{i,2};
        if n>len
            val(len+1:n) = []; % remove extra, normally for char
        elseif n<len
            val(n+1:len) = 0; % pad 0, normally for char
        end
        fwrite(fid, val, C0{i,3});
    end
    
    % Write nii ext: extension is in hdr
    for i = 1:nExt % nExt may be 0
        fwrite(fid, nii.ext(i).esize, 'int32');
        fwrite(fid, nii.ext(i).ecode, 'int32');
        fwrite(fid, nii.ext(i).edata, 'uint8');
    end
    
    if ~isNii
        fclose(fid); % done with .hdr
        fid = fopen([fname '.img'], 'w');
    end

    % Write nii image
    fwrite(fid, nii.img, fmt);
    fclose(fid); % all written

    % gzip if asked
    if do_gzip
        if isNii
            gzipOS([fname '.nii']);
        else
            gzipOS([fname '.hdr']); % better not to compress .hdr
            gzipOS([fname '.img']);
        end
    end
    
elseif strcmpi(cmd, 'hdr')
    if nargin<2, error('nii_tool(''%s'') needs second input', cmd); end
    if ~ischar(varargin{1})
        error('nii_tool(''hdr'') needs nii file name as second input'); 
    end
    
    fname = nii_name(varargin{1}, '.hdr'); % get .hdr if it is .img
    [fid, clnObj, niiVer] = fopen_nii(fname); %#ok<ASGLU>
    varargout{1} = read_hdr(fid, niiVer, C, fname);
   
elseif any(strcmpi(cmd, {'ext' 'img' 'load'})) 
    if nargin<2, error('nii_tool(''%s'') needs second input', cmd); end
    if ischar(varargin{1})
        fname = nii_name(varargin{1}, '.hdr');
    elseif isstruct(varargin{1}) && isfield(varargin{1}, 'file_name')
        fname = varargin{1}.file_name;
    else        
        error(['nii_tool(''%s'') needs a file name or hdr struct from ' ...
            'nii_tool(''hdr'') as second input'], cmd); 
    end
    
    [fid, clnObj, niiVer, isNii] = fopen_nii(fname); %#ok<ASGLU>
    nii.hdr = read_hdr(fid, niiVer, C, fname);
    
    if strcmpi(cmd, 'ext') || strcmpi(cmd, 'load') 
        if ~isempty(nii.hdr.extension) && nii.hdr.extension(1)
            nii.ext = read_ext(fid, nii.hdr);
            if strcmpi(cmd, 'ext')
                varargout{1} = nii.ext;
                return; 
            end
        elseif strcmpi(cmd, 'ext')
            varargout{1} = []; 
            return; 
        end
    end
    
    if strcmpi(cmd, 'load') || strcmpi(cmd, 'img')
        if ~isNii % close .hdr file, and open .img file
            fname = nii_name(fname, '.img');
            [fid, clnObj] = fopen_nii(fname, nii.hdr.machine); %#ok<NASGU>
        end
        nii.img = read_img(fid, nii.hdr, para);
        if strcmpi(cmd, 'img')
            varargout{1} = nii.img;
        else % load
            varargout{1} = nii;
        end
    end
elseif strcmpi(cmd, 'RGBStyle')
    styles = {'afni' '' 'mricron' 'fsl'};
    curStyle = styles{para.rgb_dim};
    if nargin<2, varargout{1} = curStyle; return; end % query only
    irgb = varargin{1};
    if isempty(irgb), irgb = 1; end % default as 'afni'
    if ischar(irgb)
        if strncmpi(irgb, 'fsl', 3), irgb = 4;
        elseif strncmpi(irgb, 'mricron', 4), irgb = 3;
        else irgb = 1;
        end
    end
    if ~any(irgb == [1 3 4])
        error('nii_tool(''RGBStyle'') can have 1, 3, or 4 as second input'); 
    end
    if nargout, varargout{1} = curStyle; end % return old one
    para.rgb_dim = irgb;
elseif strcmpi(cmd, 'cat3D')
    if nargin<2, error('nii_tool(''%s'') needs second input', cmd); end
    fnames = varargin{1};
    if ischar(fnames) % guess it is like run1*.nii
        f = dir(fnames);
        f = sort({f.name});
        fnames = strcat([fileparts(fnames) '/'], f);
    end
    
    n = length(fnames);
    if n<2 || ~iscellstr(fnames)
        error('Invalid input for nii_tool(''cat3D''): %s', varargin{1});
    end

    nii = nii_tool('load', fnames{1}); % all for first file
    nii.img(:,:,:,2:n) = 0; % pre-allocate
    % For now, omit all consistence check between files
    for i = 2:n, nii.img(:,:,:,i) = nii_tool('img', fnames{i}); end
    varargout{1} = img2datatype(nii, para); % update dim
elseif strcmpi(cmd, 'default')
    flds = {'version' 'rgb_dim'}; % may add more in the future
    for i = 1:length(flds), val.(flds{i}) = para.(flds{i}); end
    if nargin<2, varargout{1} = val; return; end % query only
    if nargout, varargout{1} = val; end % return old val
    in2 = varargin;
    if ~isstruct(in2), in2 = struct(in2{:}); end
    nam = fieldnames(in2);
    for i = 1:length(nam)
        if strcmpi(nam{i}, 'rgb_dim')
            nii_tool('RGBstyle', in2.(nam{i}));
            continue;
        end
        ind = strcmpi(nam{i}, flds);
        if isempty(ind), continue; end
        para.(flds{ind}) = in2.(nam{i});
    end
    if val.version ~= para.version, C = niiHeader(para.version); end
    
    fname = [fileparts(which(mfilename)) '/nii_tool_para.mat'];
    fid = fopen(fname, 'w'); % check writing permission
    if fid<1
        fname = [getenv('HOME') '/nii_tool_para.mat'];
    else
        fclose(fid);
    end
    try
        save(fname, '-struct', 'para', flds{:}); % overwrite
    catch me
        fprintf(2, 'Failed to save default parameters for future use.\n');
        fprintf(2, '%s\n', me.message);
    end
elseif strcmpi(cmd, 'func_handle') % make a local function avail to outside 
    varargout{1} = eval(['@' varargin{1}]);
else
    error('Invalid command for nii_tool: %s', cmd);
end
% End of main function

%% Subfunction: all nii header in the order in NIfTI-1/2 file
function [C, para] = niiHeader(niiVer)
rgb_dim = 1;
if nargin<1 || isempty(niiVer)
    fname = [getenv('HOME') '/nii_tool_para.mat'];
    if ~exist(fname, 'file')
        fname = [fileparts(which(mfilename)) '/nii_tool_para.mat'];
    end

    try
        val = load(fname);
        niiVer = val.version;
        rgb_dim = val.rgb_dim;
    catch
        niiVer = 1;
    end
end

if niiVer == 1
    C = {
    % name              len  format     value           offset
    'sizeof_hdr'        1   'int32'     348             0
    'data_type'         10  'char*1'    ''              4                                          
    'db_name'           18  'char*1'    ''              14
    'extents'           1   'int32'     16384           32
    'session_error'     1   'int16'     0               36
    'regular'           1   'char*1'    'r'             38
    'dim_info'          1   'uint8'     0               39
    'dim'               8   'int16'     ones(1,8)       40
    'intent_p1'         1   'single'    0               56
    'intent_p2'         1   'single'    0               60
    'intent_p3'         1   'single'    0               64
    'intent_code'       1   'int16'     0               68
    'datatype'          1   'int16'     0               70
    'bitpix'            1   'int16'     0               72
    'slice_start'       1   'int16'     0               74
    'pixdim'            8   'single'    zeros(1,8)      76
    'vox_offset'        1   'single'    0               108
    'scl_slope'         1   'single'    1               112
    'scl_inter'         1   'single'    0               116
    'slice_end'         1   'int16'     0               120
    'slice_code'        1   'uint8'     0               122
    'xyzt_units'        1   'uint8'     0               123
    'cal_max'           1   'single'    0               124
    'cal_min'           1   'single'    0               128
    'slice_duration'    1   'single'    0               132
    'toffset'           1   'single'    0               136
    'glmax'             1   'int32'     0               140
    'glmin'             1   'int32'     0               144
    'descrip'           80  'char*1'    ''              148
    'aux_file'          24  'char*1'    ''              228
    'qform_code'        1   'int16'     0               252
    'sform_code'        1   'int16'     0               254
    'quatern_b'         1   'single'    0               256
    'quatern_c'         1   'single'    0               260
    'quatern_d'         1   'single'    0               264
    'qoffset_x'         1   'single'    0               268
    'qoffset_y'         1   'single'    0               272
    'qoffset_z'         1   'single'    0               276
    'srow_x'            4   'single'    [1 0 0 0]       280
    'srow_y'            4   'single'    [0 1 0 0]       296
    'srow_z'            4   'single'    [0 0 1 0]       312
    'intent_name'       16  'char*1'    ''              328
    'magic'             4   'char*1'    ''              344
    'extension'         4   'uint8'     [0 0 0 0]       348
    };

elseif niiVer == 2
    C = {
    'sizeof_hdr'        1   'int32'     540             0
    'magic'             8   'char*1'    ''              4
    'datatype'          1   'int16'     0               12
    'bitpix'            1   'int16'     0               14
    'dim'               8   'int64'     ones(1,8)       16
    'intent_p1'         1   'double'    0               80
    'intent_p2'         1   'double'    0               88
    'intent_p3'         1   'double'    0               96
    'pixdim'            8   'double'    zeros(1,8)      104
    'vox_offset'        1   'int64'     0               168
    'scl_slope'         1   'double'    1               176
    'scl_inter'         1   'double'    0               184
    'cal_max'           1   'double'    0               192
    'cal_min'           1   'double'    0               200
    'slice_duration'    1   'double'    0               208
    'toffset'           1   'double'    0               216
    'slice_start'       1   'int64'     0               224
    'slice_end'         1   'int64'     0               232
    'descrip'           80  'char*1'    ''              240
    'aux_file'          24  'char*1'    ''              320
    'qform_code'        1   'int32'     0               344
    'sform_code'        1   'int32'     0               348
    'quatern_b'         1   'double'    0               352
    'quatern_c'         1   'double'    0               360
    'quatern_d'         1   'double'    0               368
    'qoffset_x'         1   'double'    0               376
    'qoffset_y'         1   'double'    0               384
    'qoffset_z'         1   'double'    0               392
    'srow_x'            4   'double'    [1 0 0 0]       400
    'srow_y'            4   'double'    [0 1 0 0]       432
    'srow_z'            4   'double'    [0 0 1 0]       464
    'slice_code'        1   'int32'     0               496
    'xyzt_units'        1   'int32'     0               500
    'intent_code'       1   'int32'     0               504
    'intent_name'       16  'char*1'    ''              508
    'dim_info'          1   'uint8'     0               524
    'unused_str'        15  'char*1'    0               525
    'extension'         4   'uint8'     [0 0 0 0]       540
    };
else
    error('Nifti version %g is not supported', niiVer);
end
if nargout<2, return; end

%   class      datatype bitpix  valpix
D = {
    'ubit1'     1       1       1 % neither mricron nor fsl support this
    'uint8'     2       8       1
    'int16'     4       16      1
    'int32'     8       32      1
    'single'    16      32      1
    'single'    32      64      2 % complex
    'double'    64      64      1
    'uint8'     128     24      3 % RGB
    'int8'      256     8       1
    'single'    511     96      3 % RGB, not in NIfTI standard?
    'uint16'    512     16      1
    'uint32'    768     32      1
    'int64'     1024    64      1
    'uint64'    1280    64      1
%     'float128'  1536    128     1 % long double, for 22nd century?
    'double'    1792    128     2 % complex
%     'float128'  2048    256     2 % long double complex
    'uint8'     2304    32      4 % RGBA
    };

para.format   =  D(:,1)';
para.datatype = [D{:,2}];
para.bitpix   = [D{:,3}];
para.valpix   = [D{:,4}];
para.rgb_dim  = rgb_dim; % dim of RGB/RGBA in NIfTI FILE
para.version  = niiVer;

%% Subfunction: update hdr based on image class for 'init' and 'save'
function [nii, fmt] = img2datatype(nii, para)
dim = size(nii.img);
ndim = numel(dim);
dim(ndim+1:7) = 1;

if ndim == 8 % RGB/RGBA data, we change img type to uint8 or single if needed
    valpix = dim(8);
    if valpix == 4 % RGBA
        typ = 'RGBA'; % error info only
        nii.img = uint8(nii.img); % NIfTI only support uint8 for RGBA
    elseif valpix == 3 % RGB, must be single or uint8
        typ = 'RGB';
        if max(nii.img(:))>1, nii.img = uint8(nii.img);
        else nii.img = single(nii.img);
        end
    else
        error('Color dimension must have length of 3 for RGB and 4 for RGBA');
    end

    dim(8) = []; % remove color-dim so numel(dim)=7 for nii.hdr
    ndim = find(dim>1, 1, 'last'); % update it
elseif isreal(nii.img)
    typ = 'real';
    valpix = 1; 
else
    typ = 'complex';
    valpix = 2; 
end

if islogical(nii.img), imgFmt = 'ubit1'; 
else imgFmt = class(nii.img);
end
ind = find(strcmp(para.format, imgFmt) & para.valpix==valpix);

if isempty(ind) % only RGB and complex can have this problem
    error('nii_tool does not support %s image of ''%s'' type', typ, imgFmt);
elseif numel(ind)>1 % unlikely
    error('Non-unique datatype found for %s image of ''%s'' type', typ, imgFmt);
end

fmt = para.format{ind};
nii.hdr.datatype = para.datatype(ind);
nii.hdr.bitpix = para.bitpix(ind);
nii.hdr.dim = [ndim dim];

if nii.hdr.sizeof_hdr == 348
    nii.hdr.glmax = round(double(max(nii.img(:)))); % we may remove these
    nii.hdr.glmin = round(double(min(nii.img(:))));
end

%% Subfunction: use pigz or system gzip if available (faster)
function gzipOS(fname)
persistent cmd; % command to run gzip
if isempty(cmd)
    cmd = check_gzip;
    if ischar(cmd)
    	cmd = [cmd ' -f ']; % overwrite if exist
    elseif islogical(cmd) && ~cmd
        fprintf(2, ['None of system pigz, gzip or Matlab gzip available. ' ...
            'Files are not compressed into gz.\n']);
    end
end

if islogical(cmd)
    if cmd, gzip(fname); delete(fname); end
    return;
end
if ispc
	[err, str] = system(['start "" /B ' cmd '"' fname '"']); % background
else
    [err, str] = system([cmd '"' fname '" &']);
end
% [err, str] = system([cmd '"' fname '"']);
if err, errorLog(['Error during compression: ' str]); end

% Deal with pigz/gzip on path or in nii_tool folder, and matlab gzip/gunzip
function cmd = check_gzip
% first, try system pigz
[err, ~] = system('pigz -V 2>&1');
if ~err, cmd = 'pigz -n'; return; end

% next, try pigz included with nii_tool
m_dir = fileparts(which(mfilename));
if ismac % pigz for mac is not included in the package
    fprintf(2, [' Please install pigz for fast compression: ' ...
        'http://rudix.org/packages/pigz.html\n']);
elseif ispc % rename back pigz for Windows. Renamed to trick Matlab Central
    try %#ok<TRYNC>
        fname = [m_dir '\pigz.win'];
        if exist(fname, 'file')
            movefile(fname, [m_dir '\pigz.exe'], 'f');
        end
        fname = [m_dir '\pthreadGC2.win'];
        if exist(fname, 'file')
            movefile(fname, [m_dir '\pthreadGC2.dll'], 'f');
        end
    end
end

cmd = fullfile(m_dir, 'pigz -n');
[err, ~] = system([cmd ' -V 2>&1']);
if ~err, return; end

% Third, try system gzip
[err, ~] = system('gzip -V 2>&1'); % gzip on system path?
if ~err, cmd = 'gzip -n'; return; end

% Lastly, try to use Matlab gzip/gunzip. Check only one, since they are paired
if isempty(which('gunzip')) || ~usejava('jvm')
    cmd = false; % none of de/compress tools available
    return;
end
    
cmd = true; % use slower matlab gzip/gunzip

%% Try to use in order of pigz, system gunzip, then matlab gunzip
function outName = gunzipOS(fname)
persistent cmd; % command to run gupzip
if isempty(cmd)
    cmd = check_gzip;
    if ischar(cmd)
    	cmd = [cmd ' -f -d ']; % overwrite if exist
    elseif islogical(cmd) && ~cmd
        error('None of system pigz, gunzip or Matlab gunzip is available');
    end
end
pth = tempdir; % unzip into temp dir

if islogical(cmd)
    if cmd, outName = gunzip(fname, pth); end
    return;
end

[pth1, outName, ext] = fileparts(fname);
outName = fullfile(pth, outName);
if ~strcmp([pth1 filesep], pth), copyfile(fname, [outName ext], 'f'); end
[err, str] = system([cmd '"' outName ext '"']); % overwrite if exist
if err, fprintf(2, 'Error during decompression:\n%s\n', str); end

%% subfunction: read hdr
function hdr = read_hdr(fid, niiVer, C, fname)
if niiVer>1, C = niiHeader(niiVer); end % C defaults for version 1
fseek(fid, 0, 'bof');
for i = 1:size(C,1)
    hdr.(C{i,1}) = fread(fid, C{i,2}, C{i,3})';
    if strcmp(C{i,3}, 'char*1')
        hdr.(C{i,1}) = deblank(char(hdr.(C{i,1})));
    end
end

hdr.version = niiVer; % for 'save', unless user asks to change
[~, ~, hdr.machine]= fopen(fid); % use it for .img file

[pth, nam, ext] = fileparts(fname); % fname may be .gz
if isempty(pth)
    [pth, nam, ext] = fileparts(which(fname)); % in current folder or path
end
pth = fullfile(getfield(what(pth), 'path')); % full path
hdr.file_name = fullfile(pth, [nam ext]); % fname with full path

%% subfunction: read ext, and decode it if known ecode
function ext = read_ext(fid, hdr)
ext = []; % to avoid error, such as no ext but hdr.extension(1) was set
fseek(fid, hdr.sizeof_hdr+4, 'bof'); % +4 skip hdr.extension
nEnd = hdr.vox_offset;
if nEnd == 0 % .hdr file
    nEnd = getfield(dir(fopen(fid)), 'bytes'); % total bytes of the .hdr file
end
i = 1; % nExt. It would be nice if hdr.extension(2) stores nExt
while ftell(fid) < nEnd
    esize = fread(fid, 1, 'int32'); % multiple of 16
    if isempty(esize) || mod(esize,16), break; end % just to be safe
    ext(i).esize = esize; %#ok<*AGROW>
    ext(i).ecode = fread(fid, 1, 'int32'); 
    ext(i).edata = fread(fid, ext(i).esize-8, '*uint8'); % -8 for esize & ecode

    % Decode edata if we know ecode
    if ext(i).ecode == 40 % Matlab: any kind of matlab variable
        nByte = typecast(ext(i).edata(1:4), 'int32'); % num of bytes of MAT data
        if strcmp(hdr.machine, 'ieee-be'), nByte = swapbytes(nByte); end
        tmp = [tempname '.mat']; % temp MAT file to save edata
        fid1 = fopen(tmp, 'w');
        fwrite(fid1, ext(i).edata(5:nByte+4)); % exclude padded zeros
        fclose(fid1);
        deleteMat = onCleanup(@() delete(tmp)); % delete temp file after done
        ext(i).edata_decoded = load(tmp); % load into struct
    elseif ext(i).ecode == 6 % plain text
        str = char(ext(i).edata');
        if isempty(strfind(str, 'dicm2nii.m'))
            ext(i).edata_decoded = deblank(str);
        else % created by dicm2nii.m
            ss = struct;
            ind = strfind(str, [';' char([0 10])]); % strsplit error in Octave
            ind = [-2 ind]; % -2+3=1: start of first para
            for j = 1:numel(ind)-1
                a = str(ind(j)+3 : ind(j+1));
                a(a==0) = []; % to be safe. strtrim wont remove null
                a = strtrim(a);
                if isempty(a), continue; end
                try 
                    eval(['ss.' a]); % put all into struct
                catch me
                    fprintf(2, '%s\n', me.message);
                    fprintf(2, 'Unrecognized text: %s\n', a);
                end
            end
            flds = fieldnames(ss); % make all vector column
            for j = 1:length(flds)
                val = ss.(flds{j});
                if isnumeric(val) && isrow(val), ss.(flds{j}) = val'; end
            end
            ext(i).edata_decoded = ss;
        end
    elseif ext(i).ecode == 4 % AFNI
        ext(i).edata_decoded = deblank(char(ext(i).edata)');
    elseif ext(i).ecode == 2 % dicom
        no_dicm_hdr = isempty(which('dicm_hdr'));
        if no_dicm_hdr && isempty(which('dicominfo')), return; end 
        tmp = [tempname '.dcm'];
        fid1 = fopen(tmp, 'w');
        fwrite(fid1, ext(i).edata);
        fclose(fid1);
        deleteDcm = onCleanup(@() delete(tmp));
        if no_dicm_hdr
            ext(i).edata_decoded = dicominfo(tmp);
        else
            ext(i).edata_decoded = dicm_hdr(tmp);
        end
    end
    i = i + 1;
end

%% subfunction: read img
function img = read_img(fid, hdr, para)
ind = para.datatype == hdr.datatype;
if ~any(ind)
    error('Datatype %g is not supported by nii_tool.', hdr.datatype);
end

dim = hdr.dim(2:8);
dim(hdr.dim(1)+1:7) = 1; % avoid some error in file
dim(dim<1) = 1;
n = prod(dim) * para.valpix(ind); % num of values
fseek(fid, hdr.vox_offset, 'bof');
img = fread(fid, n, ['*' para.format{ind}]); % * to keep original class

if any(hdr.datatype == [128 511 2304]) % RGB or RGBA
    j = para.rgb_dim;
    if isfield(hdr, 'data_type')
        a = hdr.data_type;
        if numel(a)>8 && strncmp(a, 'rgb_dim=', 8) % override user's RGB_dim
            a = str2double(a(9));
            if any(a == [1 3 4]), j = a; end % 4 is not possible
        end
    end
    dim = [dim(1:j-1) para.valpix(ind) dim(j:7)]; % length=8 now
    img = reshape(img, dim);
    img = permute(img, [1:j-1 j+1:8 j]); % put RGB(A) to dim8
elseif any(hdr.datatype == [32 1792]) % complex single/double
    img = reshape(img, [2 dim]);
    img = complex(permute(img(1,:,:,:,:,:,:,:), [2:8 1]), ... % real
                  permute(img(2,:,:,:,:,:,:,:), [2:8 1]));    % imag
else % all others: valpix=1
    if hdr.datatype==1, img = logical(img); end
    img = reshape(img, dim);
    % this will cause problem if user permutes img again, so better off
%     if isfield(hdr, 'data_type') && strcmp(hdr.data_type, 'rgb_dim=4') && ...
%           ((hdr.datatype==2 && any(dim(4) == [3 4])) || ... % uint8 RGB/RGBA
%           (hdr.datatype==16 && dim(4)==3)) % single RGB
%         img = permute(img, [1:3 5:8 4]); % put RGB(A) to dim8
%     end
end

%% Return requested fname with ext, useful for .hdr and .img files
function fname = nii_name(fname, ext)
[~, f, e] = fileparts(fname);
n = length(fname);
if strcmpi(e, '.gz')
    n = n - 3; % 3 is length('.gz')
    [~, ~, e] = fileparts(f); % .nii/.hdr/.img
end
if strcmpi(e, '.nii') || strcmpi(e, ext), return; end
if ~strcmpi(e, '.hdr') && ~strcmpi(e, '.img')
    error(['Invalid NIfTI file name: ' fname]); 
end
fname(n+(-3:0)) = ext; % if not return or error, change it
 
%% fopen NIfTI file, check endian, niiVer and nii/hdr by 'magic' for header file
% clnObj takes care of fclose(fid) and delete(fopen(fid)) if isGz 
function [fid, clnObj, niiVer, isNii] = fopen_nii(fname, endian)
if nargin<2, endian = 'ieee-le'; end % only provided for .img fopen
[fid, err] = fopen(fname, 'r', endian);
if fid<1, error([err ': ' fname]); end

n = fread(fid, 4, '*uint8')'; % sizeof_hdr or signature
isGz = isequal(n(1:2), [31 139]); % gz, tgz, tar file
fnameIn = fname; % for error msg
if isGz
    fclose(fid); % close .gz file
    fname = gunzipOS(fname); % guzipped file, return unzipped fname
    fid = fopen(fname, 'r', endian);
    n = fread(fid, 4, '*uint8')';
end

if nargout<3 % must return here for .img
    clnObj = onCleanup(@()closeFile(fid, isGz));
    return;
end

% NIfTI: isBigEndian = dim(1)<0 || dim(1)>7; 
if isequal(n, [92 1 0 0]) % LE 348, nifti 1 
    niiVer = 1;
elseif isequal(n, [0 0 1 92]) % BE 348, nifti 1
    niiVer = 1;
    fclose(fid);
    fid = fopen(fname, 'r', 'ieee-be');
elseif isequal(n, [28 2 0 0]) % LE 540, nifti 2
    niiVer = 2;
elseif isequal(n, [0 0 2 28]) % BE 540, nifti 2
    niiVer = 2;
    fclose(fid);
    fid = fopen(fname, 'r', 'ieee-be');
else
    fclose(fid);
    error('Not valid NIfTI file: %s', fnameIn);
end
clnObj = onCleanup(@()closeFile(fid, isGz)); % return it to auto-close file

if niiVer == 1, fseek(fid, 344, 'bof'); end
magic = fread(fid, 3, '*char')';
vStr = num2str(niiVer);
if strcmp(magic, ['n+' vStr])
    isNii = true;
elseif strcmp(magic, ['ni' vStr])
    isNii = false;
else % likely wrong magic. Warn user and use file ext for nii detection
    fprintf(2, 'Inconsistent sizeof_hdr and magic string for file %s\n', fnameIn);
    fprintf(2, 'sizeof_hdr: %g; magic: %s\n', typecast(n, 'int32'), magic);
    isNii = strcmpi(fname(end+(-3:0)), '.nii');
end

%% fclose and delete ungzipped file if isGz
function closeFile(fid, isGz)
if isGz % close fid then delete tmp
    fname = fopen(fid); % ungzipped file
    fclose(fid);
    delete(fname);
else % only close fid
    fclose(fid);
end

%% subfunction: get help for a command
function subFuncHelp(mfile, subcmd)
fid = fopen(which(mfile));
if fid<1, error(' %s not exists.', mfile); end
clnObj = onCleanup(@()fclose(fid));
while 1 % find first % line
    ln = strtrim(fgetl(fid));
    if feof(fid), fprintf(2, ' No help text found.\n'); return; end
    if ~isempty(ln) && ln(1) == '%', break; end
end

cr = char(10);
str = [ln(2:end) cr];
while 1
    ln = strtrim(fgetl(fid));
    if isempty(ln) || ln(1) ~= '%', break; end % first non % line
    str = [str ln(2:end) cr];
end

% detect topic line before formating the str: try each in order
topicChar = [cr ' - ']; % we rely on this for each topic: see help
str = strrep(str, [cr ' -- '], topicChar); % ' -- ' is also fine
ind = strfind(str, topicChar);
if isempty(ind), disp(str); return; end % no topicChar found. Show all help text

fakeChar = repmat(char(1), 1, numel(topicChar));
str = strrep(str, topicChar, fakeChar); % will restore later

% format for reliable syntax and paragraph detection (order is important):
cr1 = [cr ' ']; % cr with a space
chars = {[mfile '  ']   [mfile ' ']; % reduce multiple space after mfile to one
         [mfile ' (']   [mfile '(']; % remove space between mfile and (
         [mfile '( ']   [mfile '(']; % remove space after mfile(
         [cr '    ']    [cr char(9)]; % replace 4 space with tab for beauty
         cr1            cr; % remove space after cr
         };
for i = 1:size(chars, 1)
    while ~isempty(strfind(str, chars{i,1}))
        str = strrep(str, chars{i,1}, chars{i,2}); % regexprep error in Octave
    end
end
str = strrep(str, cr, cr1); % restore one space after cr

dashes = strfind(str, fakeChar); % index for ' - ' after formating
str = strrep(str, fakeChar, topicChar); % restore ' - '

prgrfs = strfind(str, [cr1 cr1]); % double blank lines
nTopic = numel(dashes);
topics = ones(1, nTopic+1); % syntax 'mfile(' before each '-' line
for i = 1:nTopic
    ind = strfind(str(1:dashes(i)), [mfile '(']); % syntax before ' - '
    if isempty(ind), continue; end % no syntax before ' - ', assume start with 1
    ind = find(prgrfs < ind(end), 1, 'last'); % previous paragraph
    if isempty(ind), continue; end
    topics(i) = prgrfs(ind) + 1; % start of this topic 
end
topics(nTopic+1) = length(str); % set last topic to the end

cmd = strrep(subcmd, '?', ''); % remove ? in case it is in subcmd
if isempty(cmd) % help for main function
    disp(str(1:topics(1))); % text before first topic
    return;
end

% find a topic with cmd syntax, and the syntax is prior to ' - '
cmd = sprintf('%s(''%s''', mfile, cmd);
for i = 1:nTopic
    ind = strfind(lower(str(topics(i):dashes(i))), lower(cmd));
    if ~isempty(ind) % found the syntax in the topic
        disp(str(topics(i):topics(i+1)));
        return;
    end
end

% if we reach here, no subcmd found in syntax
fprintf(2, ' Unknown command for %s: %s\n', mfile, subcmd);
