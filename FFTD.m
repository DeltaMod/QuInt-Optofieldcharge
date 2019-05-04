%%This function is written to handle transforms to, and application of dispersion, between time and frequency
% space. 
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
A.N    = length(x); %N length vector
A.t    = x; %Time Axis
A.t1   = A.t(1); A.t2 = A.t(end);
A.t_0  = mean(A.t);               % Centre Time Axis
A.Dt   = (A.t(end) - A.t(1));                       % Sampling interval
A.dt   = A.Dt/A.N ; 
A.w_t1 = 2*pi()/A.Dt; A.w_t2 = 2*pi()/A.dt;   % Min-Max Frequency from time graph

A.w_0      = w_0;

A.w    = x; %Frequency Axis (iff type = ftt)
A.w1 = A.w(1); A.w2 = A.w(end);
A.Dw   = A.w2 - A.w1;
A.dw   = A.Dw/A.N;
A.Dt_w = 2*pi()/A.w1; A.dt_w = 2*pi()/A.w2;   
A.t_w1 = -A.Dt_w/2; A.t_w2 = A.Dt_w/2;


fun_theta = @(theta) theta*pi();
fun_phiw = @(phi0, phi1, phi2, phi3,w,w_0) phi0                       +...  %phi_0 is Carrier Envelope Phase - Values of pi() 
                                               phi1.*((w-w_0))        +...  %phi_1 is the Group Delay 
                                               phi2.*((w-w_0).^2)/2   +...  %phi_2 is the Group Velocity Dispersion GVD
                                               phi3.*((w-w_0).^3)/6 ;       %phi_3 is the Third Order Dispersion 
A.phi = Phi; 
A.phiw = fun_phiw(A.phi(1), A.phi(2), A.phi(3), A.phi(4),A.w,A.w_0);
A.theta = fun_theta(Theta);

%% Calculate post dispersion delay, and how many units to shift the pulse by to center it %
A.tD    = A.phi(2); % The added time delay by phi_2
A.tcirc = round(A.N*(round(A.phi(2)/A.Dt)-A.phi(2)/A.Dt)); 
end
switch type

    case 'ttf'
    % E_t -> E_w
    A.ttt.Et = y;
    A.Et = y/max(y);
    E.ttf.t
    A.Etfft=fft(A.Et); %A.Etfft= flip(A.Etfft);
    A.wfft =linspace(A.w_t2,A.w_t1,length(y)); % Defining the frequency axis - it's "backwards" because the fft goes high to low freq
    A.ttf.w = A.wfft;
    A.Etfft = A.Etfft.*exp(-1i.*fun_phiw(A.phi(1),A.phi(2),A.phi(3),A.phi(4),A.wfft,A.w_0));
    A.ttf.Ew = A.Etfft;
    
    A.ttf.PHITEXT = sprintf(['\\phi_0 = ', num2str(A.phi(1),'%.4g'),'\n',...
                             '\\phi_1 = ', num2str(A.phi(2),'%.4g'),'\n',...
                             '\\phi_2 = ', num2str(A.phi(3),'%.4g'),'\n',...
                             '\\phi_3 = ', num2str(A.phi(4),'%.4g'),'\n']); 

    case 'ttt'
     % E_t -> E_w
    A.ttt.Et = y/max(y);
    A.Et = y/max(y);
    A.ttt.t = x;
    
    A.Etfft=fft(A.Et); %A.Etfft= flip(A.Etfft);
    A.wfft =linspace(A.w_t2,A.w_t1,length(y)); % Defining the frequency axis - it's "backwards" because the fft goes high to low freq
    A.ttt.w = A.wfft;
    A.Etfft = A.Etfft.*exp(-1i.*fun_phiw(A.phi(1),A.phi(2),A.phi(3),A.phi(4),A.wfft,A.w_0));
    A.ttt.Ew = A.Etfft;
    
    %E_w -> E_t
    A.Etifft =flip((ifft(A.Etfft))); %We flip it here, because otherwise it is backwards when transformed back
    Dt = 2*pi()/A.wfft(end); dt_w = 2*pi()/A.wfft(1);   
    t1 = -A.t_0; t2 = Dt-A.t_0;
    A.t_ifft = linspace(t1,t2,length(A.Et));
    A.ttt.Etdisp = A.Etifft;
    A.ttt.tdisp  = A.t_ifft;
    A.ttt.Etdispc = circshift(A.Etifft,A.tcirc);
    A.ttt.tdispc  = A.ttt.tdisp + A.tD;
    
    A.ttt.PHITEXT = sprintf(['\\phi_0 = ', num2str(A.phi(1),'%.4g'),'\n',...
                             '\\phi_1 = ', num2str(A.phi(2),'%.4g'),'\n',...
                             '\\phi_2 = ', num2str(A.phi(3),'%.4g'),'\n',...
                             '\\phi_3 = ', num2str(A.phi(4),'%.4g'),'\n']); 

    case 'ftt'
    %E_w -> E_t
    A.ftt.w  = x;
    A.ftt.Ew = y/max(y);
    A.Ew     = y/max(y); %Note: The current Ew equation includes the phi component, so you must change that to apply a different dispersion.
    
    A.Ewfft      = fftshift(ifft(A.Ew));   % Performing the fourier transform using nfft as the N=nfft number of points 
    A.ftt.Etdisp = A.Ewfft;
    A.ftt.tdisp  = linspace(A.t_w1,A.t_w2,length(A.ftt.Etdisp)); 
   
    %In this function, inputting phi is only required to keep track of
    %which was used. The code treats the ftt example has ALREADY HAVING
    %DISPERSION
    A.ftt.PHITEXT = sprintf(['\\phi_0 = ', num2str(A.phi(1),'%.4g'),'\n',...
                             '\\phi_1 = ', num2str(A.phi(2),'%.4g'),'\n',...
                             '\\phi_2 = ', num2str(A.phi(3),'%.4g'),'\n',...
                             '\\phi_3 = ', num2str(A.phi(4),'%.4g'),'\n']); 

    otherwise
        fprintf(2,'ERROR: The function FFTD did not complete, you did not enter a correct "type". Please use either ttf, ttt or ftt\n')
end
end