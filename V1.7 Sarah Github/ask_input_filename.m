function [R] = ask_input_filename(varargin)
% Demonstrate how to make a custom dialog box which returns information.
%
% Example:
%          T = GUI_36;  % T will be a string.
%
% Suggested exercise:  How would you modify this code so that the default
% answer 'Enter Some Data' CANNOT be returned?
%
%
% Author:  Matt Fig
% Date:  1/15/2010
%% Hubert Taieb 2019-07-30
% Modification to have the figure nicely around the fucci axes for the
% clustering separation
if isempty(varargin)
    varargin = cellstr('input_filename');
end
R = [];  % In case the user closes the GUI.
S.fh = figure('units','pixels',...
    'position',[500 500 300 100],...
    'menubar','none',...
    'name','Filename?',...
    'numbertitle','off',...
    'resize','off');

S.input = uicontrol('style','edit',...
    'units','pix',...
    'position',[10 60 220 30],...
    'string',varargin{1},...
    'callback',{@write_filename});

S.done = uicontrol('style','pushbutton',...
    'units','pix',...
    'position',[10 20 220 30],...
    'string','Done',...
    'callback',{@confirm});
uicontrol(S.input)  % Make the editbox active.
set(S.fh,'WindowKeyPressFcn',@keyPressCallback);
uiwait(S.fh)  % Prevent all other processes from starting until closed.

    function [] = write_filename(varargin)
        % Callback for the pushbutton.
        R = get(S.input,'string');
    end


    function [] = confirm(varargin)
        % Callback for the pushbutton.
         R = get(S.input,'string');
        close(S.fh);  % Closes the GUI, allows the new R to be returned.
    end

    function keyPressCallback(source,eventdata)
      % determine the key that was pressed
      keyPressed = eventdata.Key;
      if strcmpi(keyPressed,'return')
          % the key that was pressed was the space bar so set focus to 
          % the snap button
          uicontrol(S.done);
          % invoke the callback for the snap button
          confirm(S.done,[]);
      end

  end

end