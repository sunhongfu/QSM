
function varargout = UI1_2_R2S(varargin)
% UI1_2_R2S MATLAB code for UI1_2_R2S.fig
%      UI1_2_R2S, by itself, creates a new UI1_2_R2S or raises the existing
%      singleton*.
%
%      H = UI1_2_R2S returns the handle to a new UI1_2_R2S or the handle to
%      the existing singleton*.
%
%      UI1_2_R2S('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UI1_2_R2S.M with the given input arguments.
%
%      UI1_2_R2S('Property','Value',...) creates a new UI1_2_R2S or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UI1_2_R2S_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UI1_2_R2S_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UI1_2_R2S

% Last Modified by GUIDE v2.5 11-Mar-2020 15:16:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UI1_2_R2S_OpeningFcn, ...
                   'gui_OutputFcn',  @UI1_2_R2S_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before UI1_2_R2S is made visible.
function UI1_2_R2S_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UI1_2_R2S (see VARARGIN)

% Choose default command line output for UI1_2_R2S
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UI1_2_R2S wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UI1_2_R2S_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global path_mag;
[file_mag,path_mag] = uigetfile('*.dcm');
apath_mag = string([path_mag,file_mag]); 
set(handles.edit1,'string',path_mag);
 
 


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global path_ph;
[file_ph,path_ph] = uigetfile('*.dcm');
set(handles.edit2,'string',path_ph);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global path_out;
path_out = uigetdir();
set(handles.edit3,'string',path_out);


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in process.
function process_Callback(hObject, eventdata, handles)
% hObject    handle to process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%initialization
global path_mag;
global path_ph;
global path_out;
global readout_v;
global R_mask_v
global ph_unwrap_v;
global BFRM_v;
global options
options.readout    = 'unipolar';
options.R_mask     = 1;
options.fit_thr    = 40;
options.bet_thr    = 0.4;
options.bet_smooth = 2;
options.ph_unwrap  = 'prelude'; 
options.bkg_rm     = 'resharp';
options.t_svd      = 0.1;
options.smv_rad    = 3;
options.tik_reg    = 1e-3;
options.cgs_num    = 500;
options.lbv_peel   = 2;
options.lbv_tol    = 0.01;
options.tv_reg     = 5e-4;
options.tvdi_n     = 500;
%%
readout_v = get(handles.readout,'String');
readout_v_num = get(handles.readout,'Value');
options.readout    = readout_v{readout_v_num};

% options.readout    = readout_v;
options.R_mask     = R_mask_v;
options.bet_thr    = get(handles.bet_thr,'value');
options.bet_smooth = get(handles.bet_smooth,'value');
options.ph_unwrap  = ph_unwrap_v; 
options.bkg_rm     = BFRM_v;
options.t_svd      = get(handles.t_svd,'value');
options.smv_rad    = get(handles.smv_rad,'value');
options.tik_reg    = get(handles.tik_reg,'value');
options.cgs_num    = get(handles.cgs_num,'value');
options.lbv_peel   = get(handles.lbv_peel,'value');
options.lbv_tol    = get(handles.lbv_tol,'value');
options.tv_reg     = get(handles.tv_reg,'value');
options.tvdi_n     = get(handles.tvdi_n,'value');
qsm_r2s_prisma(path_mag, path_ph, path_out, options);


% --- Executes on selection change in readout.
function readout_Callback(hObject, eventdata, handles)
% hObject    handle to readout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global readout_v
val=get(handles.readout,'value');
switch val
    case 1
    readout_v = 'Unipolar';   
    case 2
    readout_v = 'Bipolar'; 
end
% Hints: contents = cellstr(get(hObject,'String')) returns readout contents as cell array
%        contents{get(hObject,'Value')} returns selected item from readout


% --- Executes during object creation, after setting all properties.
function readout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to readout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in R_mask.
function R_mask_Callback(hObject, eventdata, handles)
% hObject    handle to R_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global R_mask_v
val=get(handles.R_mask,'value');
switch val
    case 1
    R_mask_v = 'Yes';   
    case 2
    R_mask_v = 'NO'; 
end
% Hints: contents = cellstr(get(hObject,'String')) returns R_mask contents as cell array
%        contents{get(hObject,'Value')} returns selected item from R_mask


% --- Executes during object creation, after setting all properties.
function R_mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fit_thr_Callback(hObject, eventdata, handles)
% hObject    handle to fit_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fit_thr,'value',str2double(get(hObject,'string')));
% Hints: get(hObject,'String') returns contents of fit_thr as text
%        str2double(get(hObject,'String')) returns contents of fit_thr as a double


% --- Executes during object creation, after setting all properties.
function fit_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bet_thr_Callback(hObject, eventdata, handles)
% hObject    handle to bet_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.bet_thr,'value',str2double(get(hObject,'string')));
% Hints: get(hObject,'String') returns contents of bet_thr as text
%        str2double(get(hObject,'String')) returns contents of bet_thr as a double


% --- Executes during object creation, after setting all properties.
function bet_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bet_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bet_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to bet_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.bet_smooth,'value',str2double(get(hObject,'string')));

% Hints: get(hObject,'String') returns contents of bet_smooth as text
%        str2double(get(hObject,'String')) returns contents of bet_smooth as a double


% --- Executes during object creation, after setting all properties.
function bet_smooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bet_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ph_unwrap.
function ph_unwrap_Callback(hObject, eventdata, handles)
% hObject    handle to ph_unwrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ph_unwrap_v
val=get(handles.ph_unwrap,'value');
switch val
    case 1
    ph_unwrap_v = 'prelude';   
    case 2
    ph_unwrap_v = 'bestpath'; 
end
% Hints: contents = cellstr(get(hObject,'String')) returns ph_unwrap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ph_unwrap


% --- Executes during object creation, after setting all properties.
function ph_unwrap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ph_unwrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in default_p.
function default_p_Callback(hObject, eventdata, handles)
% hObject    handle to default_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global path_mag;
global path_ph;
global path_out;
qsm_r2s_prisma(path_mag, path_ph, path_out, []);

% --- Executes on selection change in BFRM.
function BFRM_Callback(hObject, eventdata, handles)
% hObject    handle to BFRM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BFRM_v
val=get(handles.BFRM,'value');
switch val
    case 1
    BFRM_v = 'pdf';   
    case 2
    BFRM_v = 'sharp'; 
    case 3
    BFRM_v = 'resharp';
    case 4
    BFRM_v = 'esharp';
    case 5
    BFRM_v = 'lbv';
end
% Hints: contents = cellstr(get(hObject,'String')) returns BFRM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BFRM


% --- Executes during object creation, after setting all properties.
function BFRM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BFRM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_svd_Callback(hObject, eventdata, handles)
% hObject    handle to t_svd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.t_svd,'value',str2double(get(hObject,'string')));

% Hints: get(hObject,'String') returns contents of t_svd as text
%        str2double(get(hObject,'String')) returns contents of t_svd as a double


% --- Executes during object creation, after setting all properties.
function t_svd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_svd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smv_rad_Callback(hObject, eventdata, handles)
% hObject    handle to smv_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.smv_rad,'value',str2double(get(hObject,'string')));

% Hints: get(hObject,'String') returns contents of smv_rad as text
%        str2double(get(hObject,'String')) returns contents of smv_rad as a double


% --- Executes during object creation, after setting all properties.
function smv_rad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smv_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tik_reg_Callback(hObject, eventdata, handles)
% hObject    handle to tik_reg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.tik_reg,'value',str2double(get(hObject,'string')));
% Hints: get(hObject,'String') returns contents of tik_reg as text
%        str2double(get(hObject,'String')) returns contents of tik_reg as a double


% --- Executes during object creation, after setting all properties.
function tik_reg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tik_reg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cgs_num_Callback(hObject, eventdata, handles)
% hObject    handle to cgs_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cgs_num,'value',str2double(get(hObject,'string')));
% Hints: get(hObject,'String') returns contents of cgs_num as text
%        str2double(get(hObject,'String')) returns contents of cgs_num as a double


% --- Executes during object creation, after setting all properties.
function cgs_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cgs_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lbv_peel_Callback(hObject, eventdata, handles)
% hObject    handle to lbv_peel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.lbv_peel,'value',str2double(get(hObject,'string')));
% Hints: get(hObject,'String') returns contents of lbv_peel as text
%        str2double(get(hObject,'String')) returns contents of lbv_peel as a double


% --- Executes during object creation, after setting all properties.
function lbv_peel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbv_peel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lbv_tol_Callback(hObject, eventdata, handles)
% hObject    handle to lbv_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.lbv_tol,'value',str2double(get(hObject,'string')));

% Hints: get(hObject,'String') returns contents of lbv_tol as text
%        str2double(get(hObject,'String')) returns contents of lbv_tol as a double


% --- Executes during object creation, after setting all properties.
function lbv_tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbv_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tv_reg_Callback(hObject, eventdata, handles)
% hObject    handle to tv_reg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.tv_reg,'value',str2double(get(hObject,'string')));

% Hints: get(hObject,'String') returns contents of tv_reg as text
%        str2double(get(hObject,'String')) returns contents of tv_reg as a double


% --- Executes during object creation, after setting all properties.
function tv_reg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tv_reg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tvdi_n_Callback(hObject, eventdata, handles)
% hObject    handle to tvdi_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.tvdi_n,'value',str2double(get(hObject,'string')));

% Hints: get(hObject,'String') returns contents of tvdi_n as text
%        str2double(get(hObject,'String')) returns contents of tvdi_n as a double


% --- Executes during object creation, after setting all properties.
function tvdi_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tvdi_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
