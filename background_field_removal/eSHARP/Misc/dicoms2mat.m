function mat = dicoms2mat(path_in)
%DICOMS2MAT Convert dicoms in a folder into a matrix

% MAT = DIOCMS2MAT(PATH_IN)

lists = dir( [ path_in '/*.dcm']);

for i = 1:numel(lists)
    mat(:,:,i) = double( dicomread([path_in '/' lists(i).name]) );
end
