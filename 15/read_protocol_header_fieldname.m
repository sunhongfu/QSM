function fieldvalue_to_read = read_protocol_header(filename,fieldname_to_read);
% fieldvalue_to_read = read_protocol_header(filename,fieldname_to_read,fieldtype);
%
% function to read a targetted parameter from siemens ascii protocol (mrprot.asc) header
%
% filename - name of ascii data file including the path
% fieldname_to_read - name of the field whose value is to be retreived
% fieldvalue_to_read - value of the field that is retreived
%
%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************
fieldvalue_to_read = [];
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
            fieldname = strrep(fieldname,fieldname(index2+1:index3-1),num2str(str2num(fieldname(index2+1:index3-1))+1));
        end

        if strcmp(fieldname,fieldname_to_read)
            fieldvalue_to_read = fieldvalue;
            return;
        end
        
	end
end
fclose(fid);

