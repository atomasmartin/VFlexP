function fun_close_all_dialogs
%FUN_CLOSE_ALL_DIALOGS Summary of this function goes here
%   Detailed explanation goes here
% Get the handles of all open warning dialogs
warningDialogHandles = allchild(0);

% Check if each handle is a warning dialog and close it
for i = 1:length(warningDialogHandles)
    close(warningDialogHandles(i));
    
end

end

