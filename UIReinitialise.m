%% This script is used to re-initialise the UI, in case you need to reset or initialise - Run at start/reset
addpath(genpath(fileparts(which('Variables.m'))))
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex'); set(0, 'defaultTextInterpreter','latex');
set(0,'defaultAxesFontName', 'CMU Serif'); set(0,'defaultTextFontName', 'CMU Serif');

GUIHANDLES = flip(findobj('Type','uicontrol','Style','edit'));
GUIEDITNAMES = cellstr(GUIHANDLES(1).Tag);
VAREDITNAMES = cellstr(GUIHANDLES(1).Tag(3:end));
for n = 2:length(GUIHANDLES)
GUIEDITNAMES = [GUIEDITNAMES,cellstr(GUIHANDLES(n).Tag)];
VAREDITNAMES = [VAREDITNAMES,cellstr(GUIHANDLES(n).Tag(3:end))];
end
assignin('base','GUIEDITNAMES',GUIEDITNAMES);
assignin('base','VAREDITNAMES',VAREDITNAMES);

for n = 1:length(VAREDITNAMES)
assignin('base',VAREDITNAMES{n},1)
end
if exist('LOADPRESET','var') == 0
DATAIMPORT = readtable('SessionRestore.dat','Delimiter',';');
else
DATAIMPORT = readtable(LOADPRESET,'Delimiter',';');    
clear LOADPRESET
end

%This loads both new and legacy save slots
ORIGNAMEORDER = DATAIMPORT.Properties.VariableNames(2:end); 
VAR.Name = VAREDITNAMES;   
for n = 1:length(GUIEDITNAMES)
IND = find(strcmp(ORIGNAMEORDER,VAREDITNAMES{n}));
if isempty(IND) == 1
    VAR.Mag(n)    = 1;
    VAR.Handle(n) = 1;
    VAR.Dat(n)    = 1;
elseif isempty(IND) == 0
    VAR.Mag(n)    = DATAIMPORT{2,1+IND};
    VAR.Handle(n) = DATAIMPORT{1,1+IND};
    VAR.Dat(n)    = VAR.Handle(n).*VAR.Mag(n);   
end
end
    

% VAR.Handle = DATAIMPORT{1,2:end}; VAR.Mag = DATAIMPORT{2,2:end}; VAR.Dat = VAR.Handle.*VAR.Mag; 
% VAR.Name   = DATAIMPORT.Properties.VariableNames(2:end); 
VAR.Sim    = DATAIMPORT{4,1};
VAR.Plot   = DATAIMPORT{5,1};
if height(DATAIMPORT) == 5
VAR.LLP    = cellstr('None');
else 
VAR.LLP    = DATAIMPORT{6,1};
end
 assignin('base','VAR',VAR);
for n = 1:length(VAR.Dat)
feval(@()assignin('caller',VAR.Name{n},VAR.Handle(n))); assignin('base',VAR.Name{n},VAR.Dat(n))   
end


SimulationSelect = [{'Simple Photoinduced Charge'}
                    {'Dispersion and Photoinduced Charge'}];
             
PlotSelectList{1} = [{'Gaussian Laser Pulse'        }
                    {'<a^2n+1>'                     }
                    {'Photoinduced Charge'          }
                    {'<a^2n+1> Term Contributions to Charge'}];
PlotSelectList{2} = [{'Non Fourier Et(t)'           }
                    {'E_t FFT Plot'                 }
                    {'E_t Final Applied Dispersion' }
                    {'E_t(t)^2n+1 After Dispersion' }
                    {'Post Dispersion Photoinduced Charge'}
                    {'E_w IFFT Dispersion'          }
                    {'Current Pulse Overlap at given Delta t'}
                    {'Temporal Overlap Induced Charge'}
                    {'Dual Polarised Induced Charge'}];

assignin('base','SimulationSelect',SimulationSelect);       %Store all potential dropdown menu items in base
assignin('base','SimSelect', VAR.Sim);  %Set current dropdown item from dropdown in base
assignin('base','PlotSelectList',PlotSelectList);           %Store all potential dropdown menu items in base
assignin('base','PlotSelect',VAR.Plot); %Set current dropdown item from dropdown in base
SimIndex = find(ismember(SimulationSelect, VAR.Sim));
switch SimIndex
case 1
    PlotIndex = find(ismember(PlotSelectList{1}, VAR.Plot));
    case 2
        PlotIndex = find(ismember(PlotSelectList{2}, VAR.Plot));
end

assignin('base','CONSOLEHISTORY',[{''},{''},{''},{''},{''},{''},{''},{''}])

%% Setting up Initial UI Values %%
for n = 1:length(GUIEDITNAMES)
    set(handles.(GUIEDITNAMES{n}), 'String',VAR.Handle(n))
end
% set(handles.UIf_0x,              'String', f_0x)
% set(handles.UIf_0y,              'String', f_0y)
% set(handles.UIF_0x,              'String', F_0x)
% set(handles.UIF_0y,              'String', F_0y)
% set(handles.UIAeff,              'String', Aeff)
% set(handles.UIL,                 'String', L)
% set(handles.UINcycx,             'String', Ncycx)
% set(handles.UINcycy,             'String', Ncycy)
% set(handles.UIt1,                'String', t1)
% set(handles.UIt2,                'String', t2)
% set(handles.UIN,                 'String', N)
% set(handles.UIBPF,               'String', BPF);
% set(handles.UIDelt1,             'String', Delt1)
% set(handles.UIDelt2,             'String', Delt2)
% set(handles.UIORD,               'String', ORD)

set(handles.UISimulationSelector,'Value',SimIndex);                 %Sets the value of the Simulation selector to 1
set(handles.UISimulationSelector,'String',SimulationSelect);                 %Sets the handles for all Simulation items in dropdow
set(handles.UIPlotSelect,        'Value',PlotIndex);                         %Sets the value of the Plot selector to the correct index
set(handles.UIPlotSelect,        'String',PlotSelectList{SimIndex});        %Sets the handles for all graphs for said simulation in dropdown
set(handles.UISliderPos,           'Value', 1)












%% Alternate method of data importing
% fileID = fopen('DefaultPreset.txt');                     %Get ID of file, and store in memory
% PresetContent = textscan(fileID,'%s','delimiter',';') ; %Read each line of the file (Important so we can separate table rows later)
% dl = length(PresetContent{1}(1:end-4));                 %Set Datalength
% for n = 1:4
%     for m = 1:dl/4
% DATATEST{n,m} = PresetContent{1}{(n-1)*dl/4+m}; 
%     end
% end
% for n = 5:6
%     for m = 1:2
% DATATEST{n,m} = PresetContent{1}{end-4+n-5+m}; 
%     end
% end
