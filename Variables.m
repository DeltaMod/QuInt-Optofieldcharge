%% -- Constants -- &&
eps_0 = 8.85418782 * 10^-12;    %m^-3kg^-1s^4A^2    - Permittivity of free space
Chi1  = (1:4).^2 - 1;           %No units           - First order susceptibility
Chi2  = 1;                      %Unused             - Second order susceptibility
Chi3  = 2*10^-22;               %m^2V^-2            - Third order susceptibility (This is for Silica)
F_0   = linspace(0,F_0x,5000);  %Vm^-1              - Optical Field
F2_0  = 1;                      %Unused             - 
F3_0  = 1;                      %Unused             -

%% -- Variables -- %%
c        = 299792458 ;             %m/s
lambda_0 = c/f_0x;                 %m
w_0x     = 2*pi()*f_0x;            %Rad/s              - Laser Frequency
w_0y     = 2*pi()*f_0y;
t_0      = (t1+t2)/2;
F_a      = 5.36*10^10;             % Vm^-1             - Atomic Field  

%========================================================================================%
%% == Calculating the \phi_n components of \Phi(\omega-\omega_0) - CEP, GD, GVD, TOD == %%
%========================================================================================%

%Equation describing the refractive index of fused silica for wavelengths in the range of 0.21-6.7 \mu m
%Since the equations use a wide variety of units for their calculation, we
%must convert to match what we want (fs^n/mm)

fun_n    = @(lambda) sqrt( ((0.6961663*(lambda*10^6)^2)/((lambda*10^6)^2 - 0.0684043^2)) + ...
                           ((0.4079426*(lambda*10^6)^2)/((lambda*10^6)^2 - 0.1162414^2)) + ...
                           ((0.8974794*(lambda*10^6)^2)/((lambda*10^6)^2 - 9.896161^2))  +  1); 
n_ref = fun_n(lambda_0);
%Take Derivatives of the function defining n, we need them for \Phi
syms WavLen
fun_nd1 = eval(['@(WavLen)' char(diff(fun_n(WavLen)))]);
nd1     = fun_nd1(lambda_0)*(10^-6);
fun_nd2 = eval(['@(WavLen)' char(diff(fun_nd1(WavLen)))]);
nd2     = fun_nd2(lambda_0)*(10^-6);
fun_nd3 = eval(['@(WavLen)' char(diff(fun_nd2(WavLen)))]);
nd3     = fun_nd3(lambda_0)*(10^-6);
clear WavLen
                       
%Graph of equation for the relevant range above - we don't really need this though
nind = n_ref; lambrange = linspace(0.21,6.7,1000); % \mu m

%n_ref= fun_n(lambda_0);                                                                   %units none
GD    = (1/((c/n_ref)/(1 - (lambda_0/n_ref) * nd1))) / 1000;                               %units s/mm
GVD   = (lambda_0^3)/(2 * pi * c^2) * nd2 * 1000;                                          %units s^2/mm
TOD   = -((lambda_0)/(2 * pi() * c))^2 * 1/c * (3*lambda_0^2*nd2+lambda_0^3*nd3) * 1000;   %units s^3/mm


SiO2.n   = n_ref;               % Refractive Index  - CEP replaces this later. 
SiO2.CD  = L*(nd1*10^3)*2*pi;   %                   - Cromatic dispersion
SiO2.GD  = L*GD;                % s^1               - Group Delay
SiO2.GVD = L*GVD;               % s^2               - Group Velocity dispersion
SiO2.TOD = L*TOD;               % s^3               - Third Order Dispersion


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