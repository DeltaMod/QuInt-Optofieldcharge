%% The first few variables are without description and comments from the original KhurginModelling
%% For cross referencing, please look into KhurginModelling.m instead

%% -- Modelling Optically Induced Currents -- %%

if exist('UIRUN','var') == 0
    disp('UIRUN =/= 1, running with default variables.')
SimSelect = "Dispersion and Photoinduced Charge";
addpath(genpath(fileparts(which('Variables.m'))))
%% -- Constants -- &&
fprintf('Running non-UI version -- Check Variables')
eps_0    = 8.85418782 * 10^-12;    %m^-3kg^-1s^4A^2    - Permittivity of free space
Chi1     = (1:4).^2 - 1;           %No units           - First order susceptibility
Chi2     = 1;                      %Unused             - Second order susceptibility
Chi3     = 2*10^-22;               %m^2V^-2            - Third order susceptibility (This is for Silica)
F_0      = (0:0.05:2.5).*(10^10);  %Vm^-1              - Optical Field
%% -- Variable Parameters -- %%
f_0x     = 375*10^12;             %Hz                 - Laser frequency (From our lab)                    %Hz- For Fig4, use this f_0 = 1 to get the same result as they did.
f_0y     = 375*10^12*2;            %Hz                 - Laser frequency (first Harmonic) 
c        = 299792458 ;             %m/s
lambda_0 = c/f_0x;                 % m
w_0x     = 2*pi()*f_0x;           %Rad/s              - Laser Frequency
w_0y     = 2*pi()*f_0y;            %Rad/s              - Laser Frequency
Aeff     = 2.3*10^-12;             % m^2               - Effective area 
Ncycx    = 1.5;                   % no unit           - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
Ncycy    = 1.5;                    % no unit           - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
F_a      = 5.36*10^10;             % Vm^-1             - Atomic Field
L        = 0.1;                    % mm                - Length of material dispersion
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

%n_ref= fun_n(lambda_0);                                                                  %units none
GD    = 1/((c/n_ref)/(1 - (lambda_0/n_ref) * nd1)) * 1000;                                 %units s/mm
GVD   = (lambda_0^3)/(2 * pi * c^2) * nd2 * 1000;                                          %units s^2/mm
TOD   = -((lambda_0)/(2 * pi() * c))^2 * 1/c * (3*lambda_0^2*nd2+lambda_0^3*nd3) * 1000;   %units s^3/mm


SiO2.n   = n_ref;        % Refractive Index  - CEP replaces this later. 
SiO2.GD  = L*GD;         % s^2               - Group Delay
SiO2.GVD = L*GVD;        % s^2               - Group Velocity dispersion
SiO2.TOD = L*TOD;        % s^3               - Third Order Dispersion

BPF     = 0;
%% Note that lambrange is in \mu m, so if you want to find a refractive index at a wavelength lambda, simply use fun_n(lambda/10^-6)
for n = 1:length(lambrange)
n_wav(n) = real(fun_n(lambrange(n)));
end
A_t      = 1;                      % no unit           - Amplitude
A_w      = 1;                      % no unit           - Amplitude
%% -- Defining the Electromagnetic Field Transient -- %%
T_x = (Ncycx*2*pi())/w_0x;             % seconds - Pulse Duration (downconverted)  
T_y = (Ncycx*2*pi())/w_0x;             % seconds - Pulse Duration (downconverted)  
%Documentation claims that: ???t=4ln(2), that is W*T = 4*log(2) which means we get:
W_x    = 4*log(2)./T_x;
W_y    = 4*log(2)./T_y;
t1   = 0; t2 = 15*T_x; t_0 = (t2+t1)/2;
N    = 100;
t    = linspace(t1,t2,N);             % Seconds - Time
Dt   = (t2-t1);                       % Sampling interval
dt   = Dt/N ; 
w_t1 = 2*pi()/Dt; w_t2 = 2*pi()/dt;   % Min-Max Frequency from time graph

w1   = w_0x - w_0x*0.99; w2 = w_0x + 20*w_0x; 
w    = linspace(w1,w2,N);
Dw   = w2-w1;
dw   = Dw/N;
Dt_w = 2*pi()/w1; dt_w = 2*pi()/w2;   
t_w1 = -Dt_w/2; t_w2 = Dt_w/2;

Delt1  = -100*10^(-15); Delt2 = 100*10^(-15);

Delt = linspace(Delt1,Delt2,N);
end



%% New a_fun code @(t,T,w_0x) that is more universal and can be integrated
fun_a = @(t,T,w_0x) (exp(-2*log(2).*t.^2/T^2).*cos(w_0x*t));
a = fun_a(t,T_x,w_0x);

%We define a new function handle a_fun2n1 to help us integrate it over time with the new ^2n+1
fun_a2n1 = @(t,T,n,w_0x) (exp(-2*log(2).*t.^2/T.^2).*cos(w_0x*t)).^(2*n+1);

%Integrate for relevant values of n

for n = 1:ORD
a2(n) = w_0x*trapz(t,fun_a2n1(t,T_x,n,w_0x));
end


%% NEW CODE BLOCK - Let's generate some ultrashort pulses:
%study [8] states that ultrashort pulses can be defined in: 
% Time Domain     : E(t) = E_0(t)*e^(i*(w_0x*t - psi(t)))
% Frequency Domain: E(w) = E_0(w)*e^(-i*phi(w))          (This is a fourier transform of the above)
% E_0(w) = Amplitude, phi(w) = phase

% According to
% "https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses"
% we can simplify E(t) to instead exclude the complex conjugate (c.c.)

fun_Et = @(A_t,t,t_0,T,w_0x,theta) sqrt( A_t.*exp(-log(2).*((2*(t-t_0))/(T)).^2)).*exp(-1i*(w_0x*(t-t_0)+theta)); %+ c.c.
fun_Ew = @(A_w,w,W,w_0x,psi) sqrt( A_w.*exp(-log(2).*((2*(w-w_0x))  /(W)).^2)).*exp(-1i.*psi);
                                   
% Function for phase, theta = phase shift, t = time, set to ensure the phase starts at 0                                   
%fun_theta = @(theta,t) theta*pi() - t(1)*w_0x + t*w_0x;
fun_theta = @(theta) theta*pi(); %When you find out how to define the chirp, the group velocity dispersion and other factors like the TOD, come back to this.
                                   
%Note: fun_phiw describes phi(w-w_0x), spectral phase. It has several
%components that influence where the pulse ends up.

%% Determining  SiO2.phi manually %%
%Note: k = 2*pi/lambda, and w = 2*pi()*f, f = c/lambda => w = 2*pi*c/lambda  => lambda = 2*pi*c/w   => k = 2*pi*lambda/(2*pi*c)
%Note that vg = c/(n(w)* w dn/dw) where n(w)*w dn/dw = ng.
%Differentiating the lambda version using : syms lambda; diff(symbols, lambda)
% k_lambda    =  2*pi()/lambda
% k_lambda'   = -(2*pi)/lambda^2
% k_lambda''  =  (4*pi)/lambda^3
% k_lambda''' = -(12*pi)/lambda^4
%But I am not so sure how this, in itself, would help?
% k_w on the other hand is different?:
%k_w = w/vp =  w/= ; 

fun_k = @(w,vp) w/vp;
fun_klambda = @(lambda) 2*pi()/lambda;

k(1) = w_0x/SiO2.GD;
k(2) = c/SiO2.GD;
%online sources state that k' = 1/vg where vg = c/n_g

SiO2.phi = [pi() SiO2.GD SiO2.GVD SiO2.TOD 1]; 

fun_phiw = @(phi_0, phi_1, phi_2, phi_3,w,w_0x) phi_0                 +...    %phi_0 is Carrier Envelope Phase - Values of pi() 
                                               phi_1.*((w-w_0x))      +...    %phi_1 is the Group Delay 
                                               phi_2.*((w-w_0x).^2)/2 +...    %phi_2 is the Group Velocity Dispersion GVD
                                               phi_3.*((w-w_0x).^3)/6 ;       %phi_3 is the Third Order Dispersion 
              

%We have to think for a moment about what the phase phi actually is. If we
%have an oscillation at w_0x = something, that means the phase must complete
%on the same order as the frequency.
%Note: We use T and W here, but in most documentation these are: \delta t and \delta \omega 

%% Define separate Et and Ew series based on input variables
theta = fun_theta(0);
phiw = fun_phiw(SiO2.phi(1),SiO2.phi(2),SiO2.phi(3),SiO2.phi(4),w,w_0x);

Et = fun_Et(A_t,t,t_0,T_x,w_0x,theta); 
Ew = fun_Ew(A_w,w,W_x,w_0x,phiw);
%Note: These are related to each other by: \delta w*\delta t =4ln(2), that is W*T = 4*log(2) which means we get: W  = 4*log(2)/T;
%T = (Ncycx*2*pi())/w_0x => W = 4*log(2)*w_0x/(Ncycx*2*pi()) 
PROGBAR = waitbar(0.1,'Calculating $E_t(t)$'); %pause(0.2);


%% if PlotSelect == "Non Fourier Et(t)"
% plot(t-t_0,real(Et)); xlabel('time [s]'); ylabel('E_t(t)'); fprintf(['\n','Plotting [',PlotSelect{1},']']);
% end




%%-- Time to test conversions between the two --%%
%% -- E_t -> FFT = E_w -> Band Filter + IFFT = E_t -- %% 

waitbar(0.2,PROGBAR,'Calculating fft($E_t$)');%pause(0.2);

Estr = FFTD(t,Et,w_0x,'ttt',SiO2.phi,0); %This generates:  | Estr.ttt.Et & Estr.ttt.t | Estr.ttt.Ew & Estr.ttf.w | Estr.ttt.Etdisp & Estr.ttt.tdisp|

%%


%%Bandpass Filter
if BPF > 0
waitbar(0.3,PROGBAR,'Applying Bandpass Filter')
for n = 1:length(Estr.ttt.Ew)
    if abs(Estr.ttt.Ew(n))<max(abs(Estr.ttt.Ew))/BPF
Estr.ttt.Ew(n) = 0;
    end
end
Estr2 = FFTD(Estr.ttt.w,Estr.ttt.Ew,w_0x,'ftt',SiO2.phi,0);
Estr.ttt.Etdisp = Estr2.ftt.Etdisp; Estr.ttt.tdisp = Estr2.ftt.tdisp; 
clear Estr2
end


%% if PlotSelect == "E_t FFT Plot"
%plot(Estr.ttt.w,abs(Estr.ttt.Ew)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
%end
%%

waitbar(0.3,PROGBAR,'Calculating Dispersed $E_t(t)$');

%% if PlotSelect == "E_t Final Applied Dispersion"
%plot(Estr.ttt.tdisp,real(Estr.ttt.Etdisp));grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
%end

%% -- Post Dispersion <a^2n+1> and Photoinduced Charge -- %%
%After we get back Etifft, we now have a new value for a(t), and thus need to calculate a new <a^2n+1>

%% -- Calculating <a^2n+1> for a dispersed pulse -- %%
waitbar(0.4,PROGBAR,'Calculating Dispersed Photoinduced Charge');
for n = 1:ORD
    a2disp(n) = w_0x*trapz(Estr.ttt.tdisp,real(Estr.ttt.Etdisp).^(2*n+1));
end
%% -- Equation 17 -- %%
%fun_Q2 = @(F_0x,F_a,a2disp,Aeff) eps_0*F_0x.*(F_0x/F_a).^2.*(a2disp(1) +(F_0x/F_a).^2*a2disp(2)...
%                                   +(F_0x/F_a).^4*a2disp(3)+(F_0x/F_a).^6*a2disp(4)+(F_0x/F_a).^8*a2disp(5))*Aeff;
%Q = fun_Q2(F_0,F_a,a2disp,Aeff);
Q = fun_Q(F_0,F_a,a2disp,Aeff,ORD);  


Etf2n(1,:) = real(Estr.ttt.Etdisp);
for n = 1:4
Etf2n(n+1,:) = real(Estr.ttt.Etdisp).^(2*n+1);
end
%% if PlotSelect == "E_t(t)^2n+1 After Dispersion"
%plot(Estr.ttt.tdisp,Etf2n(:,:),'LineWidth',1.3)
%end


%% if PlotSelect == "Post Dispersion Photoinduced Charge"
%plot(F_0,abs(Q));
%figure(66)
%hold on
%plot(F_0,abs(Q));
%grid on
%xlabel('Optical Field [Vm^-^1]')
%ylabel('Photinduced Charge [C]')
% set(gca, 'YScale', 'log')
% plot([0,2.5*10^10],    [1.6*10^-19 1.6*10^-19]'r--')
% legend('Ncycx = 1.0', 'Ncycx = 1.1', 'Ncycx = 1.2', 'Ncycx = 1.3', 'Ncycx = 1.4', 'Ncycx = 1.5', 'Ncycx = 1.6', 'Ncycx = 1.7', 'Ncycx = 1.8', 'Ncycx = 1.9', 'Ncycx = 2.0', 'e')
%end





%% -- E_w -> IFFT = E_t, just to show that the dispersion applied is equivalent 

%Redundant Code
%phiw  = fun_phiw(SiO2.phi(1),SiO2.phi(2),SiO2.phi(3),SiO2.phi(4),w,w_0x); % Assigning the dispersion component constants, remember that SiO2.phi(1) represents \phi_0
%Ew    = fun_Ew(A_w,w,W,w_0x,phiw);   % Using the Ew function to define the frequency space dispersed wave
%nfft  = 2^nextpow2(length(Ew));     % Determining required length (nfft) of frequency axis: Padding by 2^P, where 2^(P-1)<length(Ew) and 2^P>lenght(Ew) (if length(Ew=1000), then P = 10 for instance since 2^10 = 1024>1000
%Ewfft = fftshift(ifft(Ew,nfft));    % Performing the fourier transform using nfft as the N=nfft number of points 
%PHITEXT = sprintf(['\\phi_0 = ', num2str(SiO2.phi(1),'%.4g'),'\n',...
%                   '\\phi_1 = ', num2str(SiO2.phi(2),'%.4g'),'\n',...
%                   '\\phi_2 = ', num2str(SiO2.phi(3),'%.4g'),'\n',...
%                   '\\phi_3 = ', num2str(SiO2.phi(4),'%.4g'),'\n']);
               
%% Alternative Code %%
Eftt  = FFTD(w,Ew,w_0x,'ftt',SiO2.phi,0);

 

%% Alternate Code
waitbar(0.45,PROGBAR,'Calculating $E_w(t)$');
Estw = FFTD(w,Ew,w_0x,'ftt',SiO2.phi,0);
%% if PlotSelect ==  "E_w IFFT Dispersion"
%plot(Estw.ftt.tdisp,real(Estw.ftt.Etdisp)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
%end

%% -- Dual Pulses with Temporal Delay Delt -- %%
%To start off, the temporal delay \Delta t (We will call it Delt in the
%code) is the delay between two short pulses. Namely:
%The driving pulse (F_0 = Fx0) that we will call: Ed(t) (The oscillation direction is between plates)
%The Injection Pulse (F0y) that we will call:     Ei(t) (The oscillation direction is perpendicular to the plates)
%For all intents and purposes, these are called a(t), and were called this
%in reference in the photoinducedchargecode.
%The equation for the new vector potential <a^{2n+1}>(Delt) is now given by:

%<a^{2n+1}>(Delt) = w_0x integral(@(time) a(t)*(a(t-Delt))^2n,t(1),t(N))
%Note here, that a(t) represents the Driving pulse (F_0) and a^2n represents the Injection pulse F0y


%We start by defining the two pulses (To start off, let's make them identical)
%%fun_Et(A_t,t,T,w_0x,theta)
%%fun_Ew(A_w,w,W,w_0x,psi) - You need this->Et to do this I think...
%w_d = w_0x; w_i = w_0y; % w_d, the weaker driving pulse (w_0x or w_0xx), and w_i, the stronger injection pulse (w_0y, or w_i)
%W_d = W  ; W_i = W_y; %Respective FWHM for each




%% Dual Polarisation Equation %%
%fun_QDelt2 = @(F_0x,F_0y,a2n,Aeff) eps_0.*F_0x.*(F_0y/F_a).^2*(1/3 * a2n(1)+...
%                                                              1/5 * a2n(2).*(F_0y./F_a).^2+...
%                                                              1/7 * a2n(3).*(F_0y./F_a).^4+...
%                                                              1/9 QPol2(n) = fun_QDelt(F_0(end)/12.5,F_0(end)/1.25,F_a,a2pol(n,:),Aeff,ORD)  * a2n(4).*(F_0y./F_a).^6+...
%                                                              1/11* a2n(5).*(F_0y./F_a).^8)*Aeff; 
                                                     
anim = 1;

N_Delt = Delt2 - Delt1 + 1;
Deltn  = linspace(Delt1,Delt2,N_Delt);
E_td = FFTD(t,fun_Et(A_t,t,t_0,T_x,w_0x,0),w_0x,'ttt',SiO2.phi,0);
E_ti = FFTD(t,fun_Et(A_t,t,t_0,T_y,w_0y,0),w_0y,'ttt',SiO2.phi,0);
Delt = linspace(Delt1*dt,Delt2*dt,N_Delt);
waitbar(0.5+0.5*0/N,PROGBAR,['Integrating to determine $\left<a^2n+1\right>$','   n/N$\_$Delt = ',num2str(0),'/',num2str(N_Delt)]);


for n = 1:N_Delt
if rem(n,round(N_Delt/100,-1)) == 0
waitbar(0.5+0.5*n/N_Delt,PROGBAR,['Integrating to determine $\left<a^2n+1\right>$','   n/N$\_$Delt = ',num2str(n),'/',num2str(N_Delt)]);
%if N/n = A*N/100 1000/200
end
E_td.ttt.Etdshift = circshift(E_td.ttt.Etdisp,Deltn(n)); E_td.ttt.tdshift = circshift(E_td.ttt.tdisp,Deltn(n)); 
%E_wd = fun_Ew(A_w,w ,W,w_0x,fun_phiw(SiO2.phi(1),SiO2.phi(2)+Delt(n),SiO2.phi(3),SiO2.phi(4),w,w_0x));
%E_wi = fun_Ew(A_w,w,W_y,w_0y,fun_phiw(SiO2.phi(1),SiO2.phi(2),         SiO2.phi(3),SiO2.phi(4),w,w_0y));
%E_td = FFTD(w,E_wd,'ftt',[SiO2.phi(1),SiO2.phi(2)+n*1e-15,SiO2.phi(3),SiO2.phi(4)],0);
%E_ti = FFTD(w,E_wi,'ftt',SiO2.phi,0);
Edires = E_td.ttt.Etdshift+E_ti.ttt.Etdisp;

for m = 1:ORD
    a2harm(n,m) = w_0x*trapz(E_td.ttt.tdisp,real(Edires).^(2*m+1));
    a2pol(n,m)  = w_0x*trapz(E_td.ttt.tdisp,real(E_td.ttt.Etdshift).*real(E_ti.ttt.Etdisp).^(2*m)); 
end


QHarm(n) = fun_Q(F_0(end),F_a,a2harm(n,:),Aeff,ORD);
QPol(n)  = fun_QDelt(F_0(end)/12.5,F_0(end)/1.25,F_a,a2pol(n,:),Aeff,ORD);
 
%% In case you want animation of the transformation here :)


if anim == 0
    figure(67); clf;
subplot(3,1,1)
plot(E_ti.ttt.tdisp,real(E_ti.ttt.Etdisp))
ylim([-1.3 1.3])
hold on
plot(E_td.ttt.tdisp,real(E_td.ttt.Etdshift))
ylim([-1.3 1.3])
subplot(3,1,2)
plot(E_td.ttt.tdisp,real(Edires(n,:)))
ylim([-1.3 1.3])
subplot(3,1,3)
plot(Delt(1:n),QHarm(1:n))
drawnow

figure(68)
plot(real(E_td.ttt.Etdshift))
drawnow
end
end


%% if PlotSelect == "Temporal Overlap Induced Charge"
%plot(Delt,real(QHarm(:,:))); fprintf(['\n','Plotting [',PlotSelect{1},']']);
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])
%end

%% if PlotSelect == "Dual Polarised Induced Charge"
%plot(Delt,real(QPol(:,:)));fprintf(['\n','Plotting [',PlotSelect{1},']']);
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])
%end

%Old code used to plot all of this IN-SCRIPT. Now we plot it in a separate
%script to make it easier to change/add things (without having to do
%everything twice

waitbar(1,PROGBAR,'DONE!'); %pause(0.5);
delete(PROGBAR)
run UIGraphPlotter

THESISGRAPHS = 0;
if THESISGRAPHS == 1 %This is for when you wish to reproduce graphs for your paper, keep these constant so you don't need to type them in again
close all;
%% Plotting different phi_n dispersion graphs, we do this to show exactly how the different components interact
FIG = figure(1); clf; FIG.Name = 'E_w->E_t For Different \phi_n'; FIG.Position = FigPos(2,2);
     %Note, the last number helps 
for n = 1
    DISP.A = [0 0 0 0 0]; 
    DISP.B = [30 0 0 0 1];
    DISP.C = [0 20*(10^-15) 0 0 2];
    DISP.D = [0 0 25*(10^-15)^2 0 3];
    DISP.E = [0 0 0 60*(10^-15)^3 4]; 
    SiO2.phi = (n)*DISP.E; %For testing purposes
    SiO2.phi(5) = DISP.E(5);
phiw  = fun_phiw(SiO2.phi(1),SiO2.phi(2),SiO2.phi(3),SiO2.phi(4),w,w_0x); % Assigning the dispersion component constants, remember that SiO2.phi(1) represents \phi_0
Ew    = fun_Ew(A_w,w,W_x,w_0x,phiw);   % Using the Ew function to define the frequency space dispersed wave
nfft  = 2^nextpow2(length(Ew));     % Determining required length (nfft) of frequency axis: Padding by 2^P, where 2^(P-1)<length(Ew) and 2^P>lenght(Ew) (if length(Ew=1000), then P = 10 for instance since 2^10 = 1024>1000
Ewfft = fftshift(ifft(Ew,nfft));    % Performing the fourier transform using nfft as the N=nfft number of points 

PHITEXT = sprintf(['\\phi_0 = ', num2str(SiO2.phi(1),'%.4g'),'\n',...
                   '\\phi_1 = ', num2str(SiO2.phi(2),'%.4g'),'\n',...
                   '\\phi_2 = ', num2str(SiO2.phi(3),'%.4g'),'\n',...
                   '\\phi_3 = ', num2str(SiO2.phi(4),'%.4g'),'\n']); 
               cla reset;
plot(tfft,real(Ewfft)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
legend(PHITEXT,'Location','northwest')
xlabel('Time [s]')
ylabel('Amplitude')
xlim(10^-13*[-1 1])    
end    
end