function AssignBase(hObject,handles,ctrl)
switch ctrl
    
    case 1
    
    [file,path_mag] = uigetfile('*.dcm');
    path_mag = file_path;
    return
    case 2
    [file,file_path] = uigetfile('*.dcm');
    path_ph = file_path;
    return
end