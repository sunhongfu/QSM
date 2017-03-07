function isFieldResult = myisfield (inStruct, fieldName)
%MYISFIELD
% myisfield( Obj, fieldName)
% Obj is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
% -> Returns TRUE if fieldName exists
% -> Returns FALSE otherwise
%
% from http://www.mathworks.com/matlabcentral/fileexchange/36862-viewprofiles/content/viewProfiles/myIsField.m
%
% Copyright (c) 2012, Brandon Armstrong
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:

%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
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

