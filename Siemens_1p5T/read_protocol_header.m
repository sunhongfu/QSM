function header=read_protocol_header(filename);
% function header=read_protocol_header(filename);
%
% function to read parameters from siemens ascii protocol header
% (e.g., MrProt.asc or meas.asc for VA15 and VA21, respectively)
% and output the information in a matlab structure variable "header"

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

fid=fopen(filename);
while 1
    tline = fgetl(fid);
    if ~ischar(tline); fclose(fid); return; end
    index=findstr(tline,'=');
    if ~isempty(index)
        fieldname=deblank(tline(1:findstr(tline,'=')-1));
        fieldvalue=deblank(tline(findstr(tline,'=')+1:end));
        if max(isletter(fieldvalue))==0
            fieldvalue=str2num(fieldvalue);
        end
        fieldname=strrep(fieldname,'[','{');
        fieldname=strrep(fieldname,']','}');
        index2=findstr(fieldname,'{');
        index3=findstr(fieldname,'}');
        if ~isempty(index2)
            fieldname = strrep(fieldname,fieldname(index2+1:index3-1),...
                num2str(str2num(fieldname(index2+1:index3-1))+1));
        end
        eval(['header.',fieldname,'= fieldvalue;']);
	end
end
fclose(fid);

