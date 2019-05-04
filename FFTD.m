%%This function is written to handle transforms to, and application of 
% dispersion between time and frequency domain.
%  
%
%----------------------------------USAGE----------------------------------%
%                     A = FFTD(x,y,w_0,'type',Phi,Theta)
% x is your x axis (time or frequency axis, select a proper type to match)
%
% y is your intensity/amplitude 
%
% 'type' is a string that indicates which transform you want to do:
%                        |  ttf = time to frequency   |
%                        |  ttt = time to time        |
%                        |  ftt = frequency to time   |
%
%---------------------------External Variables----------------------------%
% Theta, is a single value phase shift in the time domain. Set to 0 usually
%
% Phi is a 4 length vector containing phi_0, phi_1, phi_2, phi_3. You can
% define them as (x,y[phi_0,phi_1,phi_2,phi_3])
%
% phi_0 changes pulse phase, phi_1 changes pulse time "location", 
% phi_2 changes pulse width, phi_3 alters pulse shape in time 
%
%
%
%------------------------------Output Format------------------------------%
%
% All variables will be in A, as per: [A] = FFTD(x,y,w_0,'type',Phi,Theta)
% the actual output will be given as:
% A.type.variable - e.g. A.ttt.Et, gives the non dispersed time axis in the
% type = ttt, time to time dispersion. For the dispersion, you must write
% A.ttt.Etdisp. The same logic goes for all variables:
%  -----------------------------------------------------------------------%
%  | A.ttf.Et & A.ttf.t | A.ttf.Ew & A.ttf.w |                           | 
%  |--------------------|--------------------|---------------------------|
%  | A.ttt.Et & A.ttt.t | A.ttt.Ew & A.ttt.w | A.ttt.Etdisp & A.ttt.tdisp|
%  |--------------------|--------------------|---------------------------|
%  |                    | A.ftt.Ew & A.ftt.w | A.ftt.Etdisp & A.ftt.tdisp|
%  -----------------------------------------------------------------------
% Note: using ftt assumes you've already applied dispersion. Inputting Phi
% is just to print PHITEXT, and keep track of which dispersion you used!
function A = FFTD(x,y,w_0,type,Phi,Theta)
%% -- Variables -- %%
if issorted(real(y)) == 1
    fprintf(2,['ERROR: Your "y" input is a sorted vector. ',...
               'Are you sure you`ve not accidentally given your time/frequency axis? \n',...
               '       If this is not the case, change issorted(y) to == 0 \n'])
else
if mean(strcmp(type,{'ttt';'ttf';'ftt'})) > 0
%Temporal Variables
A.N    = length(x);                                 % N length vector
A.t    = x;                                         % Time Axis
A.t1   = A.t(1); A.t2 = A.t(end);                   % Start and end time
A.t_0  = mean(A.t);                                 % Centre Time Axis (mean of axis)
A.Dt   = (A.t(end) - A.t(1));                       % Sampling interval
A.dt   = A.Dt/A.N ;                                 % Smallest Time interval
A.w_t1 = 2*pi()/A.Dt; A.w_t2 = 2*pi()/A.dt;         % Min-Max Frequency from time graph

%Frequency Variables
A.w    = x;                                         % Frequency Axis (iff type = ftt)
A.w_0  = w_0;                                       % Centre Frequency
A.w1 = A.w(1); A.w2 = A.w(end);                     % Highest/Lowest Freq
A.Dw   = A.w2 - A.w1;                               % Frequency Difference
A.dw   = A.Dw/A.N;                                  % Smallest Frequency inverval
A.Dt_w = 2*pi()/A.w1; A.dt_w = 2*pi()/A.w2;         % Determining Max/min time interval from frequency
A.t_w1 = -A.Dt_w/2; A.t_w2 = A.Dt_w/2;              % Determining time axis start/end variables


fun_theta = @(theta) theta*pi();                                        % This sets phase offset
fun_phiw = @(phi0, phi1, phi2, phi3,w,w_0) phi0                       +...  %phi_0 is Carrier Envelope Phase - Values of pi() 
                                               phi1.*((w-w_0))        +...  %phi_1 is the Group Delay 
                                               phi2.*((w-w_0).^2)/2   +...  %phi_2 is the Group Velocity Dispersion GVD
                                               phi3.*((w-w_0).^3)/6 ;       %phi_3 is the Third Order Dispersion 
A.phi = Phi;                                                            % Get \phi values in a struct 
A.phiw = fun_phiw(A.phi(1), A.phi(2), A.phi(3), A.phi(4),A.w,A.w_0);    % Get dispersion factor 
A.theta = fun_theta(Theta);                                             % Get phase offset

%% Calculate post dispersion delay, and how many units to shift the pulse by to center it %
A.tD    = A.phi(2); % The added time delay by phi_2
A.tcirc = round(A.N*(round(A.phi(2)/A.Dt)-A.phi(2)/A.Dt)); 
end
switch type

    case 'ttf'
    %% ----- E_t -> E_w ----- %%
    %Getting Axis Variables
    A.ttf.Et = y;                       % Get E(t)         
    A.Et = y/max(y);                    % Normalise E(t) to a(t)
    A.ttf.t = x;                        % Get time axis
    
    % Perform fft and add dispersion
    A.Etfft=fft(A.Et);                   % fft with no padding
    A.wfft =linspace(A.w_t2,A.w_t1,A.N); % Defining the frequency axis - it's "backwards" because the fft goes high to low freq
    A.ttf.w = A.wfft;                    % Get frequency axis
    A.Etfft = A.Etfft.*exp(-1i.*fun_phiw(A.phi(1),A.phi(2),A.phi(3),A.phi(4),A.wfft,A.w_0)); %Apply dispersion
    A.ttf.Ew = A.Etfft;
    % Save text of which dispersion factor was used, useful for legend entries
    A.ttf.PHITEXT = sprintf(['\\phi_0 = ', num2str(A.phi(1),'%.4g'),'\n',...
                             '\\phi_1 = ', num2str(A.phi(2),'%.4g'),'\n',...
                             '\\phi_2 = ', num2str(A.phi(3),'%.4g'),'\n',...
                             '\\phi_3 = ', num2str(A.phi(4),'%.4g'),'\n']); 

    case 'ttt'
    %% ----- E_t -> E_w ----- %%
    %Getting Axis Variables
    A.ttt.Et = y/max(y);                % Get E(t)
    A.Et = y/max(y);                    % Normalise E(t) to a(t)
    A.ttt.t = x;                        % Get time axis
    
    % Perform fft and add dispersion
    A.Etfft=fft(A.Et);                   % fft with no padding
    A.wfft =linspace(A.w_t2,A.w_t1,A.N); % Defining the frequency axis - it's "backwards" because the fft goes high to low freq
    A.ttt.w = A.wfft;                    % Get frequency axis  
    A.Etfft = A.Etfft.*exp(-1i.*fun_phiw(A.phi(1),A.phi(2),A.phi(3),A.phi(4),A.wfft,A.w_0)); %Apply dispersion
    A.ttt.Ew = A.Etfft;                  % Get E_w(t) = fft(E_t(t))
    
    %% ---- E_w -> E_t ---- %%
    A.Etifft =flip((ifft(A.Etfft)));    %We flip Ew before ifft, because otherwise it is backwards when transformed back
    Dt = 2*pi()/A.wfft(end); dt_w = 2*pi()/A.wfft(1); %Extract max/min time    
    t1 = -A.t_0; t2 = Dt-A.t_0;         % Set t1/t2
    A.t_ifft = linspace(t1,t2,A.N);     % Get time axis
    A.ttt.Etdisp = A.Etifft;            % Get dispersed Et
    A.ttt.tdisp  = A.t_ifft;            % Get new time axis (ideal case, is the same as input)
    A.ttt.Etdispc = circshift(A.Etifft,A.tcirc); % Central pulse
    A.ttt.tdispc  = A.ttt.tdisp + A.tD;          % Time shifted axis - for centralised pulse
    % Save text of which dispersion factor was used, useful for legend entries
    A.ttt.PHITEXT = sprintf(['\\phi_0 = ', num2str(A.phi(1),'%.4g'),'\n',...
                             '\\phi_1 = ', num2str(A.phi(2),'%.4g'),'\n',...
                             '\\phi_2 = ', num2str(A.phi(3),'%.4g'),'\n',...
                             '\\phi_3 = ', num2str(A.phi(4),'%.4g'),'\n']); 

    case 'ftt'
    %% ---- E_w -> E_t ---- %%
    %Getting Axis Variables
    A.ftt.w  = x;        % Get frequency axis
    A.ftt.Ew = y/max(y); % Get normalised frequency domain
    A.Ew     = y/max(y); %Note: The current Ew equation includes the phi component, so you must change that to apply a different dispersion.
    
    % Perform ifft 
    A.Ewfft      = fftshift(ifft(A.Ew));        % since we start with frequency, we need fftshift
    A.ftt.Etdisp = A.Ewfft;                     % Get dispersed Et
    A.ftt.tdisp  = linspace(A.t_w1,A.t_w2,A.N); % Get time axis
   
    %In this function, inputting phi is only required to keep track of
    %which was used. The code treats the ftt example has ALREADY HAVING
    %DISPERSION - You can modify this code easily for yourself, but this
    %was done to keep it simple!
    A.ftt.PHITEXT = sprintf(['\\phi_0 = ', num2str(A.phi(1),'%.4g'),'\n',...
                             '\\phi_1 = ', num2str(A.phi(2),'%.4g'),'\n',...
                             '\\phi_2 = ', num2str(A.phi(3),'%.4g'),'\n',...
                             '\\phi_3 = ', num2str(A.phi(4),'%.4g'),'\n']); 

    otherwise
        fprintf(2,'ERROR: The function FFTD did not complete, you did not enter a correct "type". Please use either ttf, ttt or ftt\n')
end
end