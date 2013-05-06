function isFieldResult = myisfield (inStruct, fieldName)
%MYISFIELD
% myisfield( Obj, fieldName)
% Obj is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
% -> Returns TRUE if fieldName exists
% -> Returns FALSE otherwise
isFieldResult = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
if(strcmp(f{i},strtrim(fieldName)))
isFieldResult = 1;
return;
elseif isstruct(inStruct(1).(f{i}))
isFieldResult = myisfield(inStruct(1).(f{i}), fieldName);
if isFieldResult
return;
end
end
end

