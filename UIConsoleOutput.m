%% This script handles the UI callouts - As to reduce clutter

CONSOLEHISTORY = evalin('base','CONSOLEHISTORY'); CONSOLEHISTORY{length(CONSOLEHISTORY)+1}  = CONSOLECALLBACK; assignin('base','CONSOLEHISTORY',CONSOLEHISTORY); 
PADDING = [{'    '},{'    '},{'    '},{'    '},{'    '},{'    '},{'    '},{'>>'}];
CMDOUTUI = strcat(PADDING,CONSOLEHISTORY(end-7:end));
set(handles.UIConsoleOutput,'String',CMDOUTUI)
fprintf(['\n',CONSOLECALLBACK,'\n']);