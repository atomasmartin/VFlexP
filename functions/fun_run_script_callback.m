function fun_run_script_callback(gcb, file)

if ~bdIsLibrary(bdroot)
    if nargin == 1
        [file, path] = uigetfile('*.m','Select file');
        full_filename = [path file];
    else
        full_filename = file;
    end

    if file ~= 0
        try
            [path,filename,ext] = fileparts(full_filename);
            % Clear all tables
            old_path = pwd;
                cd(path)
                evalin("base",filename);
                mdlWks = get_param(bdroot(gcb), 'ModelWorkspace');
                evalin(mdlWks, filename);
                cd(old_path)
                
        catch ME
            str_name = string({ME.stack.name}');
            str_line = string({ME.stack.line}');
            errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
            errordlg(char(errortext));
        end
    end
end
end