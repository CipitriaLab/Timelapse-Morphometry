function out = isgnu()
% out = isgnu()
% More robust way to know if computer is gnu or not 
syst = computer();
if contains(syst, 'WIN')
    out = false; %REAL CASE 
%     out = true; % When i am on window and want the gnu path
elseif contains(syst, 'GLN') || contains(syst, 'MAC')
    out = true;
%     out = false; % When i am on window and want the gnu path
end


% path_test = pwd;
% 
% if strcmp(path_test(1:5),'/usr/') || strcmp(path_test(1:5),'/run/')
%     out = true;
% else 
%     out = false;
% end


end
