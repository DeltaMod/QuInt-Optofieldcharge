%% -- Modelling Optically Induced Currents -- %%
%Note: This code runs "Simplified Photoinduced charge, and does not take
%into account dispersion
    
if exist('UIRUN','var') == 0
addpath(genpath(fileparts(which('Variables.m'))))
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex'); set(0, 'defaultTextInterpreter','latex');
set(0,'defaultAxesFontName', 'CMU Serif'); set(0,'defaultTextFontName', 'CMU Serif');

disp('UIRUN =/= 1, running with default variables.')
    addpath(genpath(fileparts(which('Variables.m'))))
 SimSelect = "Simple Photoinduced Charge";
   % close all;
%% -- Constants -- &&
eps_0 = 8.85418782 * 10^-12;            %m^-3kg^-1s^4A^2    - Permittivity of free space
F_0  = linspace(0,2.5,1000).*(10^10);  %Vm^-1              - Optical Field

%% -- Variable Parameters -- %%
f_0x = 375*10^12;                %Hz                 - Laser frequency (From our lab)
w_0x = 2*pi()*f_0x;              %Rad/s              - Laser Frequency
c        = 299792458 ;          %m/s
lambda_0 = c/f_0x;               % m
Aeff  = 2.3*10^-12;             %m^2                - Effective area 
Ncycx = 1.0;                    %no unit            - Number of optical cycles in the pulse?s FWHM (Khurgin uses 1.7 in the paper)
F_a   = 5.36*10^10;             %Vm^-1              - Atomic Field
%Note: On the topic of F_a. Khurgin's paper has a typo where he defines 
% F_a ~ 5.36*10^-10, where it should be F_a ~ 5.36*10^10

%% -- Defining the Electromagnetic Field Transient -- %%
%Note: Ncyc = (T*w_0x)/(2*pi()) => T = (Ncyc*2*pi())/w_0x
 
T_x   = (Ncycx*2*pi())/w_0x;        %seconds - Pulse Duration (downconverted)                 
t   = linspace(0,T_x*4,1000);       %seconds - Time
t_0 = (t(1) + t(end))/2;            %this centres pulse around zero.
N   = length(t);                    %N series reference
ORD = 5;                            %Order of terms n considered    
dt  = (t(end)-t(1))/N;              %smallest dt, scaled to N
 
end

OCA = (t-t_0)/T_x; % Optical Cycle Axis around t_0, used to make plotting easier!

%% New a_fun code @(t,T,w_0x) that is more universal and can be integrated
a_fun = @(t,t_0,T,w_0x) (exp(-2*log(2).*(OCA).^2).*cos(w_0x*(t-t_0)));
a = a_fun(t,t_0,T_x,w_0x);


%% -- Finding <a^2n+1> -- %%
%Before we start attempting this, we first have to have a look what defines <a^(2n+1)>.
%The study states: <a^(2n+1)> = w_0x*int^{t_max}_{t_min} = a^(2n+1)(t)dt

PROGBAR = waitbar(0.3,'Calculating $\left<a^2n+1\right>$');

%We define a new function handle a_fun2n1 to help us integrate it over time with the new ^2n+1
a_fun2n1 = @(t,t_0,T,n,w_0x) (exp(-2*log(2).*(t-t_0).^2/T.^2).*cos(w_0x*(t-t_0))).^(2*n+1);

%Then, we write a loop that integrates this for the different values of n
%-we use trapz instead of int(@(t),t,A,'Waypoints',t) because it's faster
%and generates the same answer.
for n = 1:ORD
a2(n) = w_0x*trapz(t,a_fun2n1(t,t_0,T_x,n,w_0x)); %Integrating to find <a^2n+1> for n = 1:ORD
end
%% -- Verifying <a^2n+1> by plotting against Ncyc for all values of n = 1,2,3,4,5 -- %%
Ncyc2 = linspace(1,3.5,2*length(t));
T2 = (Ncyc2*2*pi())./w_0x;
waitbar(0.6,PROGBAR,['Calculating $<a^2n+1>$ for $N_{cyc} \in [', num2str(Ncyc2(1)),',',num2str(Ncyc2(end)),']$']);
for n = 1:ORD
    for Nc = 1:length(T2)  
        ta2 = linspace(0,20*T2(Nc),1000); %time axis 2, because t2 interfers with the UI code
        t2_0 = mean(ta2);
        Vpota(n,Nc) = w_0x*trapz(ta2,a_fun2n1(ta2,t2_0,T2(Nc),n,w_0x));
    end
end


waitbar(0.9,PROGBAR,['Calculating Resultant Photoinduced Charge']);
%% -- Equation 17 -- %%
%Q = eps_0*F_0.*(F_0/F_a).^2.*(a2(1) +(F_0/F_a).^2*a2(2)...
                                    %+(F_0/F_a).^4*a2(3)+(F_0/F_a).^6*a2(4)+(F_0/F_a).^8*a2(5))*Aeff;
Q = fun_Q(F_0,F_a,a2,Aeff,ORD);     %Calculate Q using fun_Q

%% Find where which <a^2n+1> term becomes dominant %%

for i = 1:length(ATermsQ.(avecFN{1}))
for n = 1:length(avecFN)
TermValCurr(n) = ATermsQ.(avecFN{n})(i); %Get charge contribution from each term <a^2n+1> at F_0(i)
end
[~,Ind(i)]  = max(TermValCurr);          %Get index of which term <a^2n+1> at F_0(i) contributes most to Q
end
clear TermValCurr; 

%Function for getting values for only the greatest term contribution
for n = 1:length(avecFN)
Temp = find(Ind==n);
if isempty(Temp) == 0
    if 1<Temp(1)
    Temp = [Temp(1)-1 Temp];
    end
    ADOMN.(avecFN{n}) = Temp;
else
    ADOMN.(avecFN{n}) = [ADOMN.(avecFN{n-1})(end)]; %This ADOMN stores only the greatest contribution to Q for all F_0
end
end

a2n(1,:) = a;
for n = 1:ORD
a2n(n+1,:) = a_fun2n1(t,t_0,T_x,n,w_0x);
end

waitbar(1,PROGBAR,['DONE!']);
delete(PROGBAR)
evalin('base','UIGraphPlotter')

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

THESISGRAPHS = 0; %This is for when you wish to reproduce graphs for your paper, keep these constant so you don't need to type them in again
if THESISGRAPHS == 1
 
%% Showing the different envelopes of the laser pulse - Intensity and field envelopes to be precise
%close all;
figure(5)

    hold on; grid on;
    Plot.A(1) = plot(OCA,cos(w_0x*(t-t_0)),':','LineWidth',1.2);
    Plot.A(2) = plot(OCA,exp(-2*log(2).*(OCA).^2),'--','LineWidth',1.2); 
    Plot.A(3) = plot(OCA,a_fun(t,t_0,T_x,w_0x),'-.','LineWidth',1.2);
    Plot.A(4) = plot(OCA,exp(-4*log(2).*(OCA).^2),'LineWidth',1.2);
    Plot.A(5) = plot([-1/2,1/2],[0.5,0.5],'--');
    plot([0,0],[-1,1],'k');plot([-3,3],[0,0],'k');
    xlim([-2 2])
    xlabel('Optical Cycles $\left[\left(\frac{t-t_0}{T_{period}}\right)\right]$','FontSize',15)
    ylabel('Normalised Amplitude','FontSize',13)
    legend('Oscillation','Amplitude','$E(t)$','Intensity','$T$ (FWHM)')
grid on
set(gca,'fontsize', 15)
set(gca,'Box','on')

LEFT = 0.15; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    %Showing the requirements of each optical field
    axis tight

end
