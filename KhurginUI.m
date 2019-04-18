
function varargout = KhurginUI(varargin)
%KHURGINUI MATLAB code file for KhurginUI.fig
%      KHURGINUI, by itself, creates a new KHURGINUI or raises the existing
%      singleton*.
%
%      H = KHURGINUI returns the handle to a new KHURGINUI or the handle to
%      the existing singleton*.
%
%      KHURGINUI('Property','Value',...) creates a new KHURGINUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to KhurginUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      KHURGINUI('CALLBACK') and KHURGINUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in KHURGINUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KhurginUI

% Last Modified by GUIDE v2.5 19-Feb-2019 16:15:11

%% How to add new variable
%REPLACE = str2num(get(hObject,'String'))*10^mag;
%assignin('base','REPLACE',REPLACE);
%VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'REPLACE'));
%VAR.Handle(Index) = REPLACE/VAR.Mag(Index); VAR.Dat(Index) = REPLACE; assignin('base','VAR',VAR)
%guihandles.Var.REPLACE = REPLACE;
%CONSOLECALLBACK = ['Setting REPLACE = ', num2str(guihandles.Var.REPLACE), ' unit']; 
%run UIConsoleOutput; drawnow;


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KhurginUI_OpeningFcn, ...
                   'gui_OutputFcn',  @KhurginUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before KhurginUI is made visible.
function KhurginUI_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
evalin( 'base', 'clearvars *' )
%Ensure that UIRUN = true only when "run" is pressed
UIRUN = 0;     assignin('base','UIRUN',UIRUN); 
%Allow only for quick plotting once the simulation is finished
assignin('base','QuickGraph',0)
%% NOTE - Newer versions of Matlab support the use of double quotes("") to define strings, but this old version needs to use the string(''), so please take care!
run UIReinitialise
%% We load our variables from a file called "Default Variables", but I will rename this "Last Session" later
% The current format is : 
%Name: f_0,  f_0y, Aeff,  Ncyc, t1,    t2,    N,    BPF
%Variable Value
%Magnitude

% Choose default command line output for KhurginUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KhurginUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = KhurginUI_OutputFcn(~, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in UISimulationSelector.
function UISimulationSelector_Callback(hObject, eventdata, handles)
% hObject    handle to UISimulationSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SimSelectList = get(handles.UISimulationSelector,'String');
SimSelectVal = get(handles.UISimulationSelector,'Value');
SimSelectStr = get(hObject,'String');
SimSelect    = string(SimSelectStr{SimSelectVal});
assignin('base','SimSelect',SimSelect);
VAR = evalin('base','VAR'); VAR.Sim = {SimSelect}; assignin('base','VAR',VAR);
handles.SimSelect = SimSelect;
CONSOLECALLBACK = ['The ', SimSelect{1}, ' simulation will be run']; 
run UIConsoleOutput; drawnow;


assignin('base','QuickGraph',0);
%% We run code that changes the dropdown menu depending on which simulation we want to run %%
PSV = evalin('base','PlotSelectList');
if SimSelectVal == 1
set(handles.UIPlotSelect,'Value',1);             %Resets value to top item in dropdown
assignin('base','PlotSelect',string(PSV{1}{1})); %Set current dropdown item from dropdown in base
set(handles.UIPlotSelect,'String',PSV{1});       %Sets the handles for all graphs for said simulation in dropdown

elseif SimSelectVal == 2
set(handles.UIPlotSelect,'Value',1);             %Resets value to top item in dropdown
assignin('base','PlotSelect',string(PSV{2}{1})); %Set current dropdown item from dropdown in base
set(handles.UIPlotSelect,'String',PSV{2});       %Sets the handles for all graphs for said simulation in dropdown

end
% Hints: contents = cellstr(get(hObject,'String')) returns UISimulationSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item  from UISimulationSelector


% --- Executes during object creation, after setting all properties.
function UISimulationSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UISimulationSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in UIPlotSelect.
function UIPlotSelect_Callback(hObject, eventdata, handles)
% hObject    handle to UIPlotSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PlotSelectList = get(handles.UIPlotSelect,'String');
PlotSelectVal = get(handles.UIPlotSelect,'Value');
PlotSelectStr = get(hObject,'String');
PlotSelect    = string(PlotSelectStr{PlotSelectVal});
assignin('base','PlotSelect',PlotSelect);
VAR = evalin('base','VAR'); VAR.Plot = {PlotSelect}; assignin('base','VAR',VAR);
handles.PlotSelect = PlotSelect;
CONSOLECALLBACK = ['The ', PlotSelect{1}, ' graph will be plotted']; 

run UIConsoleOutput; drawnow;


if evalin('base','QuickGraph') == 1
    evalin('base','UIGraphPlotter')
end
% Hints: contents = cellstr(get(hObject,'String')) returns UIPlotSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UIPlotSelect


% --- Executes during object creation, after setting all properties.
function UIPlotSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIPlotSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UIAeff_Callback(hObject, eventdata, handles)
% hObject    handle to UIAeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UIAeff as text
%        str2double(get(hObject,'String')) returns contents of UIAeff as a double
Aeff = str2num(get(hObject,'String'))*10^-12;
assignin('base','Aeff',Aeff);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'Aeff'));
VAR.Handle(Index) = Aeff/VAR.Mag(Index); VAR.Dat(Index) = Aeff; assignin('base','VAR',VAR)
guihandles.Var.Aeff = Aeff;
CONSOLECALLBACK = ['Setting Aeff = ', num2str(guihandles.Var.Aeff), ' m^2']; 
run UIConsoleOutput; drawnow;



% --- Executes during object creation, after setting all properties.
function UIAeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIAeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function UIf_0x_Callback(hObject, eventdata, handles)
% hObject    handle to UIf_0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UIf_0x as text
%        str2double(get(hObject,'String')) returns contents of UIf_0x as a double
f_0x = str2num(get(hObject,'String'))*10^12;
assignin('base','f_0x',f_0x);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'f_0x'));
VAR.Handle(Index) = f_0x/VAR.Mag(Index); VAR.Dat(Index) = f_0x; assignin('base','VAR',VAR)
guihandles.Var.f_0x = f_0x;
% Do this to retrieve variables from your figure's workspace.
CONSOLECALLBACK = ['Setting f_0x = ', num2str(guihandles.Var.f_0x), ' s']; 
run UIConsoleOutput; drawnow;



% --- Executes during object creation, after setting all properties.
function UIf_0x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIf_0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UIf_0y_Callback(hObject, eventdata, handles)
% hObject    handle to UIf_0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UIf_0y as text
%        str2double(get(hObject,'String')) returns contents of UIf_0y as a double
f_0y = str2num(get(hObject,'String'))*10^12;
assignin('base','f_0y',f_0y);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'f_0y'));
VAR.Handle(Index) = f_0y/VAR.Mag(Index); VAR.Dat(Index) = f_0y; assignin('base','VAR',VAR)
guihandles.Var.f_0y = f_0y;
% Do this to retrieve variables from your figure's workspace.
CONSOLECALLBACK = ['Setting f_0y = ', num2str(guihandles.Var.f_0y), ' s']; 
run UIConsoleOutput; drawnow;

% --- Executes during object creation, after setting all properties.
function UIf_0y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIf_0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UINcycx_Callback(hObject, eventdata, handles)
% hObject    handle to UINcycx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UINcycx as text
%        str2double(get(hObject,'String')) returns contents of UINcycx as a double
Ncycx = str2num(get(hObject,'String'));
assignin('base','Ncycx',Ncycx);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'Ncycx'));
VAR.Handle(Index) = Ncycx/VAR.Mag(Index); VAR.Dat(Index) = Ncycx; assignin('base','VAR',VAR)
guihandles.Var.Ncycx = Ncycx;
CONSOLECALLBACK = ['Setting Ncycx = ', num2str(guihandles.Var.Ncycx), ' ']; 
run UIConsoleOutput; drawnow;



% --- Executes during object creation, after setting all properties.
function UINcycx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UINcycx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UINcycy_Callback(hObject, eventdata, handles)
% hObject    handle to UINcycy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UINcycy as text
%        str2double(get(hObject,'String')) returns contents of UINcycy as a double
Ncycy = str2num(get(hObject,'String'));
assignin('base','Ncycy',Ncycy);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'Ncycy'));
VAR.Handle(Index) = Ncycy/VAR.Mag(Index); VAR.Dat(Index) = Ncycy; assignin('base','VAR',VAR)
guihandles.Var.Ncycy = Ncycy;
CONSOLECALLBACK = ['Setting Ncycy = ', num2str(guihandles.Var.Ncycy), ' ']; 
run UIConsoleOutput; drawnow;


% --- Executes during object creation, after setting all properties.
function UINcycy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UINcycy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function UIt1_Callback(hObject, eventdata, handles)
% hObject    handle to UIt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t1 = str2num(get(hObject,'String'))*10^-15;
assignin('base','t1',t1);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 't1'));
VAR.Handle(Index) = t1/VAR.Mag(Index); VAR.Dat(Index) = t1; assignin('base','VAR',VAR)
guihandles.Var.t1 = t1;
CONSOLECALLBACK = ['Setting t1 = ', num2str(guihandles.Var.t1), ' s']; 
run UIConsoleOutput; drawnow;



% Hints: get(hObject,'String') returns contents of UIt1 as text
%        str2double(get(hObject,'String')) returns contents of UIt1 as a double


% --- Executes during object creation, after setting all properties.
function UIt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UIt2_Callback(hObject, eventdata, handles)
% hObject    handle to UIt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t2 = str2num(get(hObject,'String'))*10^-15;
assignin('base','t2',t2);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 't2'));
VAR.Handle(Index) = t2/VAR.Mag(Index); VAR.Dat(Index) = t2; assignin('base','VAR',VAR)
guihandles.Var.t2 = t2;
CONSOLECALLBACK = ['Setting t2 = ', num2str(guihandles.Var.t2), ' s']; 
run UIConsoleOutput; drawnow;



% --- Executes during object creation, after setting all properties.
function UIt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%
function UIN_Callback(hObject, eventdata, handles)
% hObject    handle to UIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nold = evalin('base','N'); assignin('base','Nold',Nold);
N = str2num(get(hObject,'String'));
assignin('base','N',N);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'N'));
VAR.Handle(Index) = N/VAR.Mag(Index); VAR.Dat(Index) = N; assignin('base','VAR',VAR)
guihandles.Var.N = N;
CONSOLECALLBACK = ['Setting N = ', num2str(guihandles.Var.N), ' ']; 
run UIConsoleOutput; drawnow;
if Nold>N
    CONSOLECALLBACK = ['Clearing Variables > ', num2str(Nold), ' ']; 
    run UIConsoleOutput; drawnow;
    assignin('base','CleanWhat','N is greater than Nold')
    evalin('base','UINCleaner')
elseif Nold<N
    CONSOLECALLBACK = ['Clearing Variables in range [N = ', num2str(N), ', Nold = ',num2str(Nold),']']; 
    run UIConsoleOutput; drawnow;
    assignin('base','CleanWhat','N is smaller than Nold')
    evalin('base','UINCleaner')
end

% --- Executes during object creation, after setting all properties.
function UIN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function UIBandPassFilter_Callback(hObject, eventdata, handles)
% hObject    handle to UIBandPassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BPF = str2num(get(hObject,'String'));
assignin('base','BPF',BPF);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'BPF'));
VAR.Handle(Index) = BPF/VAR.Mag(Index); VAR.Dat(Index) = BPF; assignin('base','VAR',VAR)
guihandles.Var.BPF = BPF;
CONSOLECALLBACK = ['Setting BPF = ', num2str(guihandles.Var.BPF), ' ']; 
run UIConsoleOutput; drawnow;



% --- Executes during object creation, after setting all properties.
function UIBandPassFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIBandPassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%
function UIDelt1_Callback(hObject, eventdata, handles)
% hObject    handle to UIDelt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UIDelt1 as text
%        str2double(get(hObject,'String')) returns contents of UIDelt1 as a double
Delt1 = round(str2num(get(hObject,'String')));
Delt2 = evalin('base','Delt2'); N = evalin('base','N');

if Delt2 - Delt1 >N
    CONSOLECALLBACK = ['Range spans more than one N, note that duplicate results will be plotted']; 
    run UIConsoleOutput; drawnow;
elseif Delt1>Delt2
    CONSOLECALLBACK = ['Delt1>Delt2 - Swapping values to allow code to run']; run UIConsoleOutput; drawnow;
    CONSOLECALLBACK = ['Delt1 = ',num2str(Delt2),' and Delt2 = ',num2str(Delt1)]; run UIConsoleOutput; drawnow;
    run UIConsoleOutput; drawnow;
    temp = Delt1; Delt1 = Delt2; Delt2 = temp;
end
assignin('base','Delt1',Delt1); assignin('base','Delt2',Delt2);
VAR = evalin('base','VAR'); Index1 = find(ismember(VAR.Name, 'Delt1')); VAR = evalin('base','VAR'); Index2 = find(ismember(VAR.Name, 'Delt2'));
VAR.Handle(Index1) = Delt1/VAR.Mag(Index1); VAR.Dat(Index1) = Delt1; VAR.Handle(Index2) = Delt2/VAR.Mag(Index2); VAR.Dat(Index2) = Delt2; 
assignin('base','VAR',VAR);
guihandles.Var.Delt1 = Delt1;
guihandles.Var.Delt2 = Delt2;
set(handles.UIDelt1,             'String', Delt1)
set(handles.UIDelt2,             'String', Delt2)
CONSOLECALLBACK = ['Setting Delt1 = ', num2str(guihandles.Var.Delt1)]; 
run UIConsoleOutput; drawnow;

% --- Executes during object creation, after setting all properties.
function UIDelt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIDelt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Coffee.
function Coffee_Callback(hObject, eventdata, handles)
cla reset;
% hObject    handle to Coffee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imshow('Coffee.jpg')

% --- Executes on button press in UIRun.
%% -- This function is where all of the code goes! -- %%
function UIRun_Callback(hObject, eventdata, handles)
cla reset;
UIRUN = 1;     assignin('base','UIRUN',UIRUN);
% hObject    handle to UIRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Here's some code to use in the future if you need to verify if a variable exists:
%%if [exist('Aeff','var'),exist('f_0x','var'),exist('Ncycx','var')] == [1,1,1]
%%    disp('GOOD')
%%end

PlotSelect = evalin('base','PlotSelect');
SimSelect = evalin('base','SimSelect');
%% -- Constants are in Variables.m, add more as needed -- %%

%% -- Variables -- %%
f_0x     = evalin('base','f_0x');   %Hz
f_0y     = evalin('base','f_0y');  %Hz
F_0x     = evalin('base','F_0x');   %Hz
F_0y     = evalin('base','F_0y');  %Hz
w_0x     = 2*pi()*f_0x;             %Rad/s              - Laser Frequency
w_0y     = 2*pi()*f_0y;            %Rad/s              - Laser Frequency
Aeff     = evalin('base','Aeff');   %m^2                - Effective area  
L        = evalin('base','L');      %mm                - Sample Thickness 
Ncycx    = evalin('base','Ncycx');  % Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
Ncycy     = evalin('base','Ncycy');  % Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
Delt1    = evalin('base','Delt1'); % Shift of t_0 start
Delt2    = evalin('base','Delt2'); % Shift of t_0 end
BPF      = evalin('base','BPF');    % No unit           - Bandpass filter MaxAmp/Filter
N        = evalin('base','N');
T_y      = (Ncycx*2*pi())/w_0x;          %seconds - Pulse Duration (downconverted)
T_y      = (Ncycy*2*pi())/w_0y;          %seconds - Pulse Duration (downconverted)
t1       = evalin('base','t1'); t2 = evalin('base','t2');
t        = linspace(t1,t2,N);    %seconds - Time
t_0      = (t1+t2)/2; assignin('base','t_0', t_0)    

evalin('base','Variables')
CONSOLECALLBACK = ['Loading Variables...']; 
run UIConsoleOutput; drawnow; pause(0.2)
run Variables


%% -- Verifying <a^2n+1> by plotting against Ncyc for all values of n = 1,2,3,4,5 -- %%
if SimSelect ==  "Simple Photoinduced Charge"   
    CONSOLECALLBACK = ['Running PhotoInducedChargetoUI.m, please wait...']; 
    run UIConsoleOutput; drawnow; plot([0,1,1,0,0],[0,0,1,1,0]); text(0.5,0.5,'LOADING...','FontSize',40,'HorizontalAlignment','center','VerticalAlignment','middle'); axis tight; drawnow;
    
    evalin('base','PhotoInducedChargetoUI')
    assignin('base','QuickGraph',1);
    CONSOLECALLBACK = ['Successfully Ran Simulation'];
run UIConsoleOutput; drawnow;
elseif SimSelect == "Dispersion and Photoinduced Charge"   
    CONSOLECALLBACK = ['Running DispersionToUI.m']; 
    run UIConsoleOutput; drawnow; plot([0,1,1,0,0],[0,0,1,1,0]); text(0.5,0.5,'LOADING...','FontSize',40,'HorizontalAlignment','center','VerticalAlignment','middle'); axis tight; drawnow;
    
    evalin('base','DispersionToUI')
    assignin('base','QuickGraph',1);
    CONSOLECALLBACK = ['Successfully Ran Simulation'];
run UIConsoleOutput; drawnow;
end
pause(0.2)
CONSOLECALLBACK = ['You may now plot another graph from the dropdown']; 
run UIConsoleOutput; drawnow;
% --- Executes on button press in UIClearAll.
function UIClearAll_Callback(hObject, eventdata, handles)
% hObject    handle to UIClearAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
VAR = evalin('base','VAR');
FileName = ('SavedPresets/SessionRestore.dat');
LName    = [sprintf('Name ; '     ),sprintf('%s ; ',VAR.Name{:})  ]; LName   = LName(1:end-3);
LHandle  = [sprintf('Handle ; '   ),sprintf('%d ; ',VAR.Handle(:))]; LHandle = LHandle(1:end-3);
LMag     = [sprintf('Magnitude ; '),sprintf('%d ; ',VAR.Mag(:))   ]; LMag    = LMag(1:end-3);
LDat     = [sprintf('Data ; '     ),sprintf('%d ; ',VAR.Dat(:))   ]; LDat    = LDat(1:end-3);
LSim     = [                        sprintf('%s ; ',VAR.Sim{:})   ]; LSim    = LSim(1:end-3);
LPlot    = [                        sprintf('%s ; ',VAR.Plot{:})  ]; LPlot   = LPlot(1:end-3);
FULLTEXT =  sprintf([LName,'\n',LHandle,'\n',LMag,'\n',LDat,'\n',LSim,'\n',LPlot]);
dlmwrite(FileName,FULLTEXT,'delimiter','')
fprintf(['\nRestoring to Current Settings: ',FileName,'\n',' = ','\n',FULLTEXT])


CONSOLECALLBACK = ['Clearing all variables...']; run UIConsoleOutput; drawnow;

evalin( 'base', 'clearvars *' )
plot(0,0)

UIRUN = 0;     assignin('base','UIRUN',UIRUN); 
assignin('base','QuickGraph',0);
%% Define the initial values that go into the UI when starting it up %%
run UIReinitialise
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
VAR = evalin('base','VAR');
FileName = ('SavedPresets/SessionRestore.dat');
LName    = [sprintf('Name ; '     ),sprintf('%s ; ',VAR.Name{:})  ]; LName   = LName(1:end-3);
LHandle  = [sprintf('Handle ; '   ),sprintf('%d ; ',VAR.Handle(:))]; LHandle = LHandle(1:end-3);
LMag     = [sprintf('Magnitude ; '),sprintf('%d ; ',VAR.Mag(:))   ]; LMag    = LMag(1:end-3);
LDat     = [sprintf('Data ; '     ),sprintf('%d ; ',VAR.Dat(:))   ]; LDat    = LDat(1:end-3);
LSim     = [                        sprintf('%s ; ',VAR.Sim{:})   ]; LSim    = LSim(1:end-3);
LPlot    = [                        sprintf('%s ; ',VAR.Plot{:})  ]; LPlot   = LPlot(1:end-3);
FULLTEXT =  sprintf([LName,'\n',LHandle,'\n',LMag,'\n',LDat,'\n',LSim,'\n',LPlot]);
dlmwrite(FileName,FULLTEXT,'delimiter','')
fprintf(['\nSaved Preset as: ',FileName,'\n',' = ','\n',FULLTEXT])
delete(hObject);



% --- Executes on button press in UISavePreset.
function UISavePreset_Callback(hObject, eventdata, handles)
% hObject    handle to UISavePreset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
VAR = evalin('base','VAR');
[FileName,PathName] = uiputfile('SavedPresets/*.txt','Select a name for your preset!');
if FileName ~= 0
LName    = [sprintf('Name ; '     ),sprintf('%s ; ',VAR.Name{:})  ]; LName   = LName(1:end-3);
LHandle  = [sprintf('Handle ; '   ),sprintf('%d ; ',VAR.Handle(:))]; LHandle = LHandle(1:end-3);
LMag     = [sprintf('Magnitude ; '),sprintf('%d ; ',VAR.Mag(:))   ]; LMag    = LMag(1:end-3);
LDat     = [sprintf('Data ; '     ),sprintf('%d ; ',VAR.Dat(:))   ]; LDat    = LDat(1:end-3);
LSim     = [                        sprintf('%s ; ',VAR.Sim{:})   ]; LSim    = LSim(1:end-3);
LPlot    = [                        sprintf('%s ; ',VAR.Plot{:})  ]; LPlot   = LPlot(1:end-3);
FULLTEXT =  sprintf([LName,'\n',LHandle,'\n',LMag,'\n',LDat,'\n',LSim,'\n',LPlot]);
dlmwrite([PathName,FileName],FULLTEXT,'delimiter','')
fprintf(['\nSaved Preset as: ',FileName,'\n',' = ','\n',FULLTEXT])
    CONSOLECALLBACK = ['Saved Preset']; run UIConsoleOutput; drawnow;
else
    CONSOLECALLBACK = ['Cancelled Save']; run UIConsoleOutput; drawnow;
end


% --- Executes on button press in UILoadPreset.
function UILoadPreset_Callback(hObject, eventdata, handles)
LOADPRESET = uigetfile('SavedPresets/*.txt','Select a Preset');
if LOADPRESET ~= 0    
PRESETCALLBACK = LOADPRESET;
CONSOLECALLBACK = ['Loading Preset...']; run UIConsoleOutput; drawnow;
run UIReinitialise
CONSOLECALLBACK = ['Preset: ',PRESETCALLBACK,' loaded']; run UIConsoleOutput; drawnow; clear PRESETCALLBACK;
else
CONSOLECALLBACK = ['No Preset Selected - Resuming using old values']; run UIConsoleOutput; drawnow;
end
% hObject    handle to UILoadPreset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function UISliderDelt_Callback(hObject, eventdata, handles)
% hObject    handle to UISliderDelt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N = evalin('base','N');
VAR = evalin('base','VAR'); Index1 = find(ismember(VAR.Name, 'Delt1')); Index2 = find(ismember(VAR.Name, 'Delt2'));
Delt1 = VAR.Handle(Index1); Delt2 = VAR.Handle(Index2); 
set(hObject,'Min',Delt1); set(hObject,'Max',Delt2); set(hObject,'SliderStep',[1/(Delt2-Delt1), 1/(Delt2-Delt1)]); 
Slidern   =  round(get(hObject,'Value'));
t1 = evalin('base','t1'); t2 = evalin('base','t2'); dt = (t2-t1)/N;
CurrDeltT = Slidern*dt;
assignin('base','Slidern',Slidern);
set(handles.UISelectedDelay,'String',num2str(CurrDeltT/1e-15))
set(handles.UICurrentDeltn,'String',num2str(Slidern))
if evalin('base','QuickGraph') == 1
    evalin('base','UIGraphPlotter')
end
%%Note: Make the text box editable for the future! We might need that!


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function UISliderDelt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UISliderDelt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end 



function UIDelt2_Callback(hObject, eventdata, handles)
% hObject    handle to UIDelt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Delt2 = round(str2num(get(hObject,'String')));
Delt1 = evalin('base','Delt1'); N = evalin('base','N');

if Delt2 - Delt1 >N
    CONSOLECALLBACK = ['Range spans more than one N, note that duplicate results will be plotted']; 
    run UIConsoleOutput; drawnow;
elseif Delt1>Delt2
    CONSOLECALLBACK = ['Delt1>Delt2 - Swapping values to allow code to run']; run UIConsoleOutput; drawnow;
    CONSOLECALLBACK = ['Delt1 = ',num2str(Delt2),' and Delt2 = ',num2str(Delt1)]; run UIConsoleOutput; drawnow;
    run UIConsoleOutput; drawnow;
    temp = Delt1; Delt1 = Delt2; Delt2 = temp;
end
assignin('base','Delt1',Delt1); assignin('base','Delt2',Delt2);
VAR = evalin('base','VAR'); Index1 = find(ismember(VAR.Name, 'Delt1')); VAR = evalin('base','VAR'); Index2 = find(ismember(VAR.Name, 'Delt2'));
VAR.Handle(Index1) = Delt1/VAR.Mag(Index1); VAR.Dat(Index1) = Delt1; VAR.Handle(Index2) = Delt2/VAR.Mag(Index2); VAR.Dat(Index2) = Delt2; 
assignin('base','VAR',VAR);
guihandles.Var.Delt1 = Delt1;
guihandles.Var.Delt2 = Delt2;
set(handles.UIDelt1,             'String', Delt1)
set(handles.UIDelt2,             'String', Delt2)
CONSOLECALLBACK = ['Setting Delt2 = ', num2str(guihandles.Var.Delt2)]; 
run UIConsoleOutput; drawnow;


% --- Executes during object creation, after setting all properties.
function UIDelt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIDelt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UIF_0x_Callback(hObject, eventdata, handles)
% hObject    handle to UIF_0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

F_0x = str2num(get(hObject,'String'))*10^10;
assignin('base','F_0x',F_0x);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'F_0x'));
VAR.Handle(Index) = F_0x/VAR.Mag(Index); VAR.Dat(Index) = F_0x; assignin('base','VAR',VAR)
guihandles.Var.F_0x = F_0x;
CONSOLECALLBACK = ['Setting F_0x = ', num2str(guihandles.Var.F_0x), ' Vm^-^1']; 
run UIConsoleOutput; drawnow;

% --- Executes during object creation, after setting all properties.
function UIF_0x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIF_0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UIF_0y_Callback(hObject, eventdata, handles)
% hObject    handle to UIF_0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
F_0y = str2num(get(hObject,'String'))*10^10;
assignin('base','F_0y',F_0y);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'F_0y'));
VAR.Handle(Index) = F_0y/VAR.Mag(Index); VAR.Dat(Index) = F_0y; assignin('base','VAR',VAR)
guihandles.Var.F_0y = F_0y;
CONSOLECALLBACK = ['Setting F_0y = ', num2str(guihandles.Var.F_0y), ' Vm^-^1']; 
run UIConsoleOutput; drawnow;

% --- Executes during object creation, after setting all properties.
function UIF_0y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIF_0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UIL_Callback(hObject, eventdata, handles)
% hObject    handle to UIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

L = str2num(get(hObject,'String'));
assignin('base','L',L);
VAR = evalin('base','VAR'); Index = find(ismember(VAR.Name, 'L'));
VAR.Handle(Index) = L/VAR.Mag(Index); VAR.Dat(Index) = L; assignin('base','VAR',VAR)
guihandles.Var.L = L;
CONSOLECALLBACK = ['Setting L = ', num2str(guihandles.Var.L), ' mm']; 
run UIConsoleOutput; drawnow;

% --- Executes during object creation, after setting all properties.
function UIL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UIGraphpopout.
function UIGraphpopout_Callback(hObject, eventdata, handles)
% hObject    handle to UIGraphpopout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

POPOUTCHECKER = evalin( 'base', 'exist(''POPOUT'',''var'') == 1' );
switch POPOUTCHECKER
    case 0
    QuickGraph = evalin('base','QuickGraph');
    if QuickGraph == 1
    set(handles.UIGraphpopout,'String','Dock Graph');  
    POPOUT = 1; assignin('base','POPOUT',POPOUT); PlotSelect = evalin('base','PlotSelect');
    CONSOLECALLBACK = ['Plotting the graph: [', PlotSelect{1}, '] in a separate window']; 
    run UIConsoleOutput; drawnow;
    evalin('base','UIGraphPlotter')
    else
    CONSOLECALLBACK = ['Please run the simulation first']; 
    run UIConsoleOutput; drawnow;   
    end
    case 1
evalin('base','clear POPOUT')
close figure 1
set(handles.UIGraphpopout,'String','Popout Graph');
end
    



function UICurrentDeltn_Callback(hObject, eventdata, handles)
% hObject    handle to UICurrentDeltn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Slidern = round(str2num(get(hObject,'String')));
assignin('base','Slidern',Slidern);
N = evalin('base','N');
t1 = evalin('base','t1'); t2 = evalin('base','t2'); dt = (t2-t1)/N;
CurrDeltT = Slidern*dt;
set(handles.UISelectedDelay,'String',num2str(CurrDeltT/1e-15))
set(handles.UICurrentDeltn,'String',num2str(Slidern))
set(handles.UISliderDelt,'Value',Slidern)
CONSOLECALLBACK = ['Setting Slidern = ', num2str(Slidern)]; 
run UIConsoleOutput; drawnow;


% Hints: get(hObject,'String') returns contents of UICurrentDeltn as text
%        str2double(get(hObject,'String')) returns contents of UICurrentDeltn as a double


% --- Executes during object creation, after setting all properties.
function UICurrentDeltn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UICurrentDeltn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
