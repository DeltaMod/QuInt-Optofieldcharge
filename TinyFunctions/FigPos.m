%%%% --------------- Screensize and Plotting Coordinates --------------- %%
%
%- This function loads your screensize, and permits you to choose 9 preset
%  locations to evenly distribute your graphs without problems. 
%- Simply give coordinates for which 9th of the screen you want it displayed 
%  in as x and y coordinates.
%
%%% ---------------------- Locations FigPos(x,y) ---------------------- %%
%
%                     (1,3)   |   (2,3)   |   (3,3)
%                     -----------------------------
%                     (1,2)   |   (2,2)   |   (3,2)
%                     -----------------------------
%                     (1,1)   |   (2,1)   |   (3,1)
%
%%% --------------------------- How to use ---------------------------- %%
%
%- Fig(n) = figure('Position',FigPos(x,y)) 
%
%
%%% Enjoy This Fun Function, because you can't spell function without fun %%

function Screenlocale = FigPos(x,y)
scrsz = get(groot,'ScreenSize'); % Finds Screensize:
FigurePos{1,1} = [0         ,   0,  scrsz(3)/3, scrsz(4)/2.45];             %Location: Bottom Left
FigurePos{2,1} = [scrsz(3)/3,   0,  scrsz(3)/3, scrsz(4)/2.45];             %Location: Bottom Middle
FigurePos{3,1} = [2*scrsz(3)/3, 0,  scrsz(3)/3, scrsz(4)/2.45];             %Location: Bottom Right 
FigurePos{1,2} = [0         , scrsz(4)/3,  scrsz(3)/3, scrsz(4)/2.45];      %Location: Middle Left
FigurePos{2,2} = [scrsz(3)/3, scrsz(4)/3,  scrsz(3)/3, scrsz(4)/2.45];      %Location: Middle Middle  
FigurePos{3,2} = [2*scrsz(3)/3, scrsz(4)/3,  scrsz(3)/3, scrsz(4)/2.45];    %Location: Middle Right
FigurePos{1,3} = [0         , scrsz(4)/2,  scrsz(3)/3, scrsz(4)/2.45];      %Location: Top Left
FigurePos{2,3} = [scrsz(3)/3, scrsz(4)/2,  scrsz(3)/3, scrsz(4)/2.45];      %Location: Top Middle
FigurePos{3,3} = [2*scrsz(3)/3, scrsz(4)/2,  scrsz(3)/3, scrsz(4)/2.45];    %Location: Top Right
Screenlocale   = FigurePos{x,y};

