%% -- Constants -- &&
eps_0 = 8.85418782 * 10^-12;    %m^-3kg^-1s^4A^2    - Permittivity of free space
Chi1  = (1:4).^2 - 1;           %No units           - First order susceptibility
Chi2  = 1;                      %Unused             - Second order susceptibility
Chi3  = 2*10^-22;               %m^2V^-2            - Third order susceptibility (This is for Silica)
F_0   = linspace(0,F_0x,1000);  %Vm^-1              - Optical Field
F2_0  = 1;                      %Unused             - 
F3_0  = 1;                      %Unused             -

%% -- Variables -- %%
c        = 299792458 ;             %m/s
lambda_0 = f_0x/c;  % m
w_0x     = 2*pi()*f_0x;             %Rad/s              - Laser Frequency
w_0y     = 2*pi()*f_0y;
t_0      = (t1+t2)/2;
F_a      = 5.36*10^10;             % Vm^-1             - Atomic Field    
SiO2.GVD = L*36.0*(10^-15)^2;        % s^2               - Group Velocity dispersion
SiO2.TOD = L*27.0*(10^-15)^3;        % s^3               - Third Order Dispersion
SiO2.n   = L*1.4533;
SiO2.ng  = L*1.4672;
fun_n    = @(lambda) sqrt( ((0.6961663*lambda.^2)/(lambda.^2 - 0.0684043^2)) + ...
                           ((0.4079426*lambda.^2)/(lambda.^2 - 0.1162414^2)) + ...
                           ((0.8974794*lambda.^2)/(lambda.^2 - 9.896161^2))  +  1);
nind = fun_n(lambda_0); lambrange = linspace(0.2,7,1000); % \mu m

A_t      = 1;                      % no unit           - Amplitude
A_w      = 1;                      % no unit           - Amplitude


%% -- Defining the Electromagnetic Field Transient -- %%

T_x = (Ncycx*2*pi())/w_0x ;             % seconds - Pulse Duration (downconverted)  
T_y = (Ncycy*2*pi())/w_0y;
%Documentation claims that: ???t=4ln(2), that is W*T = 4*log(2) which means we get:
W_x  = 4*log(2)/T_x;
W_y   = 4*log(2)/T_y;
t    = linspace(t1,t2,N);             % Seconds - Time
Dt   = (t2-t1);                       % Sampling interval
dt   = Dt/N ; 
w_t1 = 2*pi()/Dt; w_t2 = 2*pi()/dt;   % Min-Max Frequency from time graph

w1   = w_0x - w_0x*0.99; w2 = w_0x + 6*w_0x; 
w    = linspace(w1,w2,N);
Dw   = w2-w1;
dw   = Dw/N;
Dt_w = 2*pi()/w1; dt_w = 2*pi()/w2;   
t_w1 = -Dt_w/2; t_w2 = Dt_w/2;

Slidern = 0;