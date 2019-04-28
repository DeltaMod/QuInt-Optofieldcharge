%% -- Modelling Optically Induced Currents -- %%
%Note: This code runs "Simplified Photoinduced charge, and does not take
%into account dispersion
    
if exist('UIRUN','var') == 0
    disp('UIRUN =/= 1, running with default variables.')
    addpath(genpath(fileparts(which('Variables.m'))))
 SimSelect = "Simple Photoinduced Charge";
   % close all;
%% -- Constants -- &&
eps_0 = 8.85418782 * 10^-12;            %m^-3kg^-1s^4A^2    - Permittivity of free space
Chi1  = (1:4).^2 - 1;                   %No units           - First order susceptibility
Chi2  = 1;                              %Unused             - Second order susceptibility
Chi3  = 2*10^-22;                       %m^2V^-2            - Third order susceptibility (This is for Silica)
F_0  = linspace(0,2.5,1000).*(10^10);  %Vm^-1              - Optical Field
F2_0  = 1;                              %Unused             - 
F3_0  = 1;                              %Unused             -


%% -- Variable Parameters -- %%

f_0x = 375*10^12;                %Hz                 - Laser frequency (From our lab)
w_0x = 2*pi()*f_0x;              %Rad/s              - Laser Frequency
c        = 299792458 ;          %m/s
lambda_0 = c/f_0x;               % m
Aeff  = 2.3*10^-12;             %m^2                - Effective area 
Ncycx = 1.0;                    %no unit            - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
F_a   = 5.36*10^10;             %Vm^-1              - Atomic Field
%Note: On the topic of F_a. In the study, it is specifically stated that
%F_a ~ 5.36*10^-10 (that is, that it is to the power of MINUS 10. If,
%however, a 10 is substituted then we exactly replicate the results from
%the paper... Whether this is a typo for us, or Khurgin, is to be
%determined later
%% -- Defining the Electromagnetic Field Transient -- %%

%Note: Ncyc = (T*w_0x)/(2*pi()) => T = (Ncyc*2*pi())/w_0x
 
T_x   = (Ncycx*2*pi())/w_0x;          %seconds - Pulse Duration (downconverted)                 
t   = linspace(0,T_x*4,1000);       %seconds - Time
t_0 = (t(1) + t(end))/2;
N   = length(t);                  %N series reference
dt  = (t(end)-t(1))/N;

L        = 0.1;                    % mm                - Length of material dispersion
SiO2.GVD = L*36.0*(10^-15)^2;      % mm s^2            - Group Velocity dispersion
SiO2.TOD = L*27.0*(10^-15)^3;      % mm s^3            - Third Order Dispersion
SiO2.n   = 1.4533;
SiO2.ng  = 1.4672;
fun_n    = @(lambda) sqrt( ((0.6961663*lambda.^2)/(lambda.^2 - 0.0684043^2)) + ...
                           ((0.4079426*lambda.^2)/(lambda.^2 - 0.1162414^2)) + ...
                           ((0.8974794*lambda.^2)/(lambda.^2 - 9.896161^2))  +  1);
nind = fun_n(lambda_0); lambrange = linspace(0.2,7,1000); % \mu m
%% Note that lambrange is in \mu m, so if you want to find a refractive index at a wavelength lambda, simply use fun_n(lambda/10^-6)
for n = 1:length(lambrange)
n_wav(n) = real(fun_n(lambrange(n)));
end
A_t      = 1;                      % no unit           - Amplitude
A_w      = 1;                      % no unit           - Amplitude
end

OCA = (t-t_0)/T_x; % Optical Cycle Axis around t_0, used to make plotting easier!

%% New a_fun code @(t,T,w_0x) that is more universal and can be integrated
a_fun = @(t,t_0,T,w_0x) (exp(-2*log(2).*(OCA).^2).*cos(w_0x*(t-t_0)));
a = a_fun(t,t_0,T_x,w_0x);


%% -- Finding <a^2n+1> -- %%
%Before we start attempting this, we first have to have a look what defines <a^(2n+1)>.

%The study states: <a^(2n+1)> = w_0x*int^{t_max}_{t_min} = a^(2n+1)(t)dt


PROGBAR = waitbar(0.3,'Calculating <a^2n+1>'); pause(0.2);
%We define a new function handle a_fun2n1 to help us integrate it over time with the new ^2n+1
a_fun2n1 = @(t,t_0,T,n,w_0x) (exp(-2*log(2).*(t-t_0).^2/T.^2).*cos(w_0x*(t-t_0))).^(2*n+1);

%Then, we write a loop that integrates this for the different values of n
%that we need.
QTERMS = 5; %Add this to the GUI later
for n = 1:QTERMS
%a2test(n) = w_0x*integral(@(time) a_fun2n1(time,T,n,w_0x),-inf+t_0,inf-t_0,'Waypoints',t) %This gives the EXACT same answer as below, making it redundant
a2(n) = w_0x*trapz(t,a_fun2n1(t,t_0,T_x,n,w_0x));
end
%% -- Verifying <a^2n+1> by plotting against Ncyc for all values of n = 1,2,3,4,5 -- %%
Ncyc2 = linspace(1,3.5,2*length(t));
T2 = (Ncyc2*2*pi())./w_0x;
waitbar(0.6,PROGBAR,['Calculating <a^2n+1> for Ncyc \in [', num2str(Ncyc2(1)),',',num2str(Ncyc2(end)),']']); pause(0.2);
for n = 1:QTERMS
    for Nc = 1:length(T2)  
        ta2 = linspace(0,20*T2(Nc),1000); %time axis 2, because t2 interfers with the UI code
        t2_0 = mean(ta2);
        %Vpota(n,Nc) = w_0x*integral(@(time) a_fun2n1(time,T2(Nc),n,w_0x),-inf,inf,'Waypoints',t); %Again, without WAYPOINTS you will get a very bad outcome, 
        Vpota(n,Nc) = w_0x*trapz(ta2,a_fun2n1(ta2,t2_0,T2(Nc),n,w_0x));
    end
end

%% if PlotSelect == "<a^2n+1>"
%semilogy(Ncyc2,abs(Vpota(:,:)))
%end

waitbar(0.9,PROGBAR,['Calculating Resultant Photoinduced Charge']); pause(0.2);
%% -- Equation 17 -- %%
%Q = eps_0*F_0.*(F_0/F_a).^2.*(a2(1) +(F_0/F_a).^2*a2(2)...
                                    %+(F_0/F_a).^4*a2(3)+(F_0/F_a).^6*a2(4)+(F_0/F_a).^8*a2(5))*Aeff;
Q = fun_Q(F_0,F_a,a2,Aeff,QTERMS);

%% Find where which <a^2n+1> term becomes dominant %%
OLDTERM = 1;

for i = 1:length(ATermsQ.(avecFN{1}))
for n = 1:length(avecFN)
TermValCurr(n) = ATermsQ.(avecFN{n})(i);
end
[Max,Ind(i)]  = max(TermValCurr);
end
clear Max; clear TermValCurr; 

for n = 1:length(avecFN)
Temp = find(Ind==n);
if isempty(Temp) == 0
    if 1<Temp(1)
    Temp = [Temp(1)-1 Temp];
    end
    ADOMN.(avecFN{n}) = Temp;
else
    ADOMN.(avecFN{n}) = [ADOMN.(avecFN{n-1})(end)];
end
end

a2n(1,:) = a;
for n = 1:4
a2n(n+1,:) = a_fun2n1(t,t_0,T_x,n,w_0x);
end

%% if PlotSelect == "Gaussian Laser Pulse"
%plot(OCA,a2n(:,1:N),'LineWidth',1.3)
%end

%if PlotSelect == "Photoinduced Charge"
%plot(F_0x,Q);
%end
waitbar(1,PROGBAR,['DONE!']); pause(0.5);
delete(PROGBAR)
run UIGraphPlotter

EXTRAQ = 0;
if EXTRAQ == 1
    figure(5);clf;
for Nc = 1:20
%% Visualising for more values of Ncyc

Ncycx = linspace(1,1.7,20);
T_x = (Ncycx.*2*pi())/w_0x; 
for n = 1:5
a2(n) = w_0x*integral(@(time) a_fun2n1(time,T_x(Nc),n,w_0x),t(1),t(N));
end
Q(Nc,:) = eps_0*F_0.*(F_0/F_a).^2.*(a2(1) +(F_0/F_a).^2*a2(2)...
                                    +(F_0/F_a).^4*a2(3))*Aeff;
 
plot(F_0,Q(Nc,:));
grid on
hold on
xlabel('Optical Field [Vm^-^1]')
ylabel('Photinduced Current [C]')


end
end

THESISGRAPHS = 1; %This is for when you wish to reproduce graphs for your paper, keep these constant so you don't need to type them in again
if THESISGRAPHS == 0
 
%% Showing the different envelopes of the laser pulse - Intensity and field envelopes to be precise
close all;
figure;

    hold on; grid on;
    Plot.A(1) = plot(OCA,cos(w_0x*(t-t_0)),':','LineWidth',1.2);
    Plot.A(2) = plot(OCA,exp(-2*log(2).*(OCA).^2),'--','LineWidth',1.2); 
    Plot.A(3) = plot(OCA,a_fun(t,t_0,T_x,w_0x),'-.','LineWidth',1.2);
    Plot.A(4) = plot(OCA,exp(-4*log(2).*(OCA).^2),'LineWidth',1.2);
    Plot.A(5) = plot([-1/2,1/2],[0.5,0.5],'--');
    plot([0,0],[-1,1],'k');plot([-3,3],[0,0],'k');
    xlabel('Optical Cycles [(t-t_0)/T]','FontSize',13)
    ylabel('Normalised Amplitude','FontSize',13)
    legend('Oscillation','Amplitude','E(t)','Intensity','T (FWHM)')

    %Showing the requirements of each optical field

end
