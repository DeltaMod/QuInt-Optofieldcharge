%% This UI Document is used to plot the graphs, and nothing else! %%
if exist('POPOUT','var') == 1 %Old Pos [525 270 870 540] %Use for simple photo [450 270 720 540] Orig: [550 270 920 520] %[550 270 560 520]
FIG = figure(1); FIG.Renderer = 'painters'; FIG.Position = [550 270 920 520] ; FIG.Name = PlotSelect{1}; %Old Ratio: [480 270 960 540]
end

%Generate Legend Names

for n = 1:length(avecFN)
tempstr = num2str(2*n+1);
avecLEG{n} = ['$\left<a^{',tempstr,'}\right>$']; 
atLEG{n+1} = ['$a(t)^{',tempstr,'}$'];
end
atLEG{1} = '$a(t)$';
clear tempstr1; clear tempstr2; clear tempstr3; clear tempstr4; clear finstr; clear finstr1; clear finstr2;

%% PhotoinducedChargetoUI %%

if SimSelect == "Simple Photoinduced Charge"
if exist('UIRUN','var') == 0; PlotSelect = "Gaussian Laser Pulse"; FIG = figure(1);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,1); end
if PlotSelect == "Gaussian Laser Pulse"
cla reset
plot((t-t_0)/T_x,a2n(:,1:N),'LineWidth',1.3)
ylim([-1 1])
set(gca, 'XLimSpec', 'Tight');
xlabel('Optical Cycle ($\frac{t-t_0}{T}$)')
ylabel('Vector Potential (and $a^{2n+1}$)')
legend(atLEG,'Location','northeast','FontSize',16)
end

if exist('UIRUN','var') == 0; PlotSelect = "<a^2n+1>"; FIG = figure(2);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,3); end
if PlotSelect == "<a^2n+1>"
cla reset;
semilogy(Ncyc2,Vpota(:,:))
hold on
%axis([1 2.5 10^-6 10^0])

legend(avecLEG,'Location','best' )
xlabel('Number of Optical Cycles $\left(\frac{t-t_0}{T}\right)$')
ylabel('Vector Potential momenta')
end


if exist('UIRUN','var') == 0; PlotSelect = "<a^2n+1> Term Contributions to Charge";  FIG = figure(3); clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(2,3); end
if PlotSelect == "<a^2n+1> Term Contributions to Charge"
    cla reset

for n = 1:length(avecFN)
area(F_0(ADOMN.(avecFN{n})),ATermsQ.(avecFN{n})(ADOMN.(avecFN{n})),'FaceAlpha',0.3,'LineStyle','None')
hold on
end
set(gca,'ColorOrderIndex',1)
for n = 1:length(avecFN)
semilogy(F_0,ATermsQ.(avecFN{n}),'HandleVisibility','off')
hold on
end




for n = length(avecFN):-1:1
    if length(ADOMN.(avecFN{n})) > 2
        AQMax = max(ATermsQ.(avecFN{n}));
        break
    end
end

set(gca,'YScale','log')
xlim([0 max(F_0)])
ylim([10^-20 AQMax])
clear AQMax;
LEG = legend({avecLEG{1:end-round(ORD-ORD/7)},'etc'},'Location','Best');
title(LEG,sprintf('Greatest Term\n Contribution'))
xlabel('Field Strength [Vm$^{-1}$]')
ylabel('Contribution to Q by $\left<a^{2n+1}\right>$ term [C]')

end

          

if exist('UIRUN','var') == 0; PlotSelect = "Photoinduced Charge";  FIG = figure(4); clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(2,3); end
if PlotSelect == "Photoinduced Charge"
    cla reset
plot(F_0,Q/(10^-15));

xlabel('Optical Field [Vm$^{-1}$]')
ylabel('Photoinduced Charge [fC]')
end
end



%% DispersionToUI %%
if SimSelect == "Dispersion and Photoinduced Charge"
if exist('UIRUN','var') == 0; PlotSelect = "Non Fourier Et(t)"; FIG = figure(1);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(1,3); end
if PlotSelect == "Non Fourier Et(t)"
cla reset;
%subplot(1,2,1)
plot(t-t_0,real(Et)); xlabel('time [s]'); ylabel('$E_t(t)$'); fprintf(['\n','Plotting [',PlotSelect{1},']']);


end

if exist('UIRUN','var') == 0; PlotSelect = "E_t FFT Plot"; FIG = figure(2);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(2,3); end
if PlotSelect == "E_t FFT Plot"
cla reset;
plot(Estr.ttt.w,abs(Estr.ttt.Ew)); fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('$\omega$ [rad/s]')
ylabel('Amplitude')
hold on
plot([w_0x,w_0x],[0,abs(max(Estr.ttt.Ew))])
end

if exist('UIRUN','var') == 0; PlotSelect = "E_t Final Applied Dispersion"; FIG = figure(3);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end 
if PlotSelect == "E_t Final Applied Dispersion"
cla reset;
plot(Estr.ttt.tdisp/10^-15,real(Estr.ttt.Etdispc)); fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('Time [fs]')
LEG = legend(['$\phi_2$ = ',num2str(Estr.tD), ' [fs]']);
title(LEG,'Pulse Delay')
ylabel('Amplitude')

end

if exist('UIRUN','var') == 0; PlotSelect = "E_t(t)^2n+1 After Dispersion"; FIG = figure(4);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,1); end
if PlotSelect == "E_t(t)^2n+1 After Dispersion"
cla reset
plot(Estr.ttt.tdisp,Etf2n(:,:),'LineWidth',1.3)
hold on
ylim([-1 1])
xlabel('Optical Cycle $\left(\frac{t-t_0}{T}\right)$')
ylabel('Vector Potential (and $a^{2n+1}$)')
legend('$E_t(t)$','$E_t(t)^3$','$E_t(t)^5$','$E_t(t)^7$','$E_t(t)^9$','Location','best' )
 
end

if exist('UIRUN','var') == 0; PlotSelect = "Post Dispersion Photoinduced Charge";  FIG = figure(5); clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(2,3); end
if PlotSelect == "Post Dispersion Photoinduced Charge"
    cla reset
    plot(F_0,(Q)/10^-15);

xlabel('Optical Field [Vm$^{-1}$]')
ylabel('Photinduced Charge [fC]')

THESISGRAPH = 1; % Warning: This seems to be really buggy, don't expect it
%to work outright when you try -Code is sloppy!
if exist('THESISGRAPH') == 1
    NREP = 120;
    if exist('TGRPH') == 0
        TGRPH  = 1;
    end
while TGRPH < NREP
LLeg{TGRPH} = ['L = ',num2str(L*10^3),' $\mu m$'];
Lrng(TGRPH) = L;
L = Lrng(1)+(TGRPH-1)*1e-3;
Qpc(TGRPH,:) = Q;
TGRPH  = 1+TGRPH;
if TGRPH<NREP+1
run Variables.m
run DispersionToUI.m
end  

end

if exist('PLOTDONE') == 0
ff =  figure(1);
plot(F_0,(Qpc)/10^-15);
ffl = legend(LLeg); ffl.NumColumns = 2;
xlabel('Field Strength [Vm$^{-1}$]')
ylabel('Photoinduced Charge [fC]')
ff.Position = [550 270 920 520];
grid on
set(gca,'fontsize', 15)
set(gca,'Box','on')
LEFT = 0.12; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;    
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])

figure(2)
plot(Lrng*10^3,Qpc(:,end)/10^-15)
xlabel('L [$\mu m$]')
ylabel('Photoinduced Charge [fC]')
grid on
set(gca,'fontsize', 15)
set(gca,'Box','on')
axis tight
LEFT = 0.12; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;    
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])

end
PLOTDONE = 1;
end
end

if exist('UIRUN','var') == 0; PlotSelect = "E_w IFFT Dispersion"; FIG = figure(6);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(1,1); end 
if PlotSelect ==  "E_w IFFT Dispersion"
cla reset;
plot(Estw.ftt.tdisp,real(Estw.ftt.Etdisp)); fprintf(['\n','Plotting [',PlotSelect{1},']']);
legend(['$\phi_2$ = ' num2str(SiO2.phi(3),'%.4g')])
xlabel('Time [s]')
ylabel('Amplitude')
end


if exist('UIRUN','var') == 0; PlotSelect = "Current Pulse Overlap at given Delta t"; FIG = figure(7);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end
if PlotSelect == "Current Pulse Overlap at given Delta t"
SldInd = find(ismember(Deltn,Slidern));
ETD = E_td; ETI = E_ti;
ETD.ttt.Etdshift = circshift(ETD.ttt.Etdisp,Deltn(SldInd)); ETD.ttt.tdshift = circshift(ETD.ttt.tdisp,Deltn(SldInd)); 
Ediresdt = ETD.ttt.Etdshift+ETI.ttt.Etdisp;

for m = 1:ORD
    a2harmdt(m) = w_0x*trapz(ETD.ttt.tdisp,real(Ediresdt).^(2*m+1));
    a2poldt(m)  = w_0x*trapz(ETD.ttt.tdisp,real(ETD.ttt.Etdshift).*real(ETI.ttt.Etdisp).^(2*m)); 
end
QHarmdt(:) = fun_Q(F_0(end),F_a,a2harmdt(:),Aeff,ORD);
QPoldt(:)  = fun_QDelt(F_0x,F_0y,F_a,a2poldt(:),Aeff,ORD);
cla reset;

hold on
plt1 = area(ETI.ttt.tdisp,real(Ediresdt)-1,-1,'FaceColor','b'); alpha(plt1,0.5);
plt2 = plot(ETD.ttt.tdisp,real(ETD.ttt.Etdshift)-1,'r','LineWidth',1);
plt3 = plot(ETI.ttt.tdisp,real(ETI.ttt.Etdisp)-1,'g','LineWidth',1);
plt4 = plot(Delt./(N_Delt/N),QHarm/max(QHarm)+1,'m'); 
plt5 = plot(Delt(SldInd)./(N_Delt/N),QHarmdt/max(QHarm)+1,'m*');
plt6 = plot(Delt./(N_Delt/N),QPol/max(QPol)+1,'k'); 
plt7 = plot(Delt(SldInd)./(N_Delt/N),QPoldt/max(QPol)+1,'k*'); 
%plot(ETI.ttt.tdisp,real(Ediresdt),'b:','LineWidth',1.3);
fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('$\Delta t$ [s]')
ylabel('Sum of Vector Potential')
axis([ETI.ttt.tdisp(1) ETI.ttt.tdisp(end) -3 3])
hold off
end


if exist('UIRUN','var') == 0; PlotSelect = "Temporal Overlap Induced Charge"; FIG = figure(8);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end
if PlotSelect == "Temporal Overlap Induced Charge"
cla reset;
plot(Delt/10^-15,real(QHarm));fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('$\Delta t$ [fs]')
ylabel('$Q(\Delta t$ [C])')
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])


%figure(1)
%plot(1,1)
end

if exist('UIRUN','var') == 0; PlotSelect = "Pulse Temporal Overlap"; FIG = figure(9);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end
if PlotSelect == "Pulse Temporal Overlap"
cla reset;
plot(E_td.ttt.tdisp,real(E_td.ttt.Etdshift(1,:)),'-.r')
hold on
plot(E_ti.ttt.tdisp,real(E_ti.ttt.Etdisp),'--b')
ylim([-2 2])
plot(E_td.ttt.tdisp,real(Edires(1,:)),':k','LineWidth',1.2)
fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('$\Delta t$ [s]')
ylabel('Vector Potential Momenta')
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])
%figure(1)
%plot(1,1)
end

if exist('UIRUN','var') == 0; PlotSelect = "Dual Polarised Induced Charge"; FIG = figure(10);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end
if PlotSelect == "Dual Polarised Induced Charge"
cla reset;
plot(Delt/10^-15,real(QPol)/10^-15);fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('$\Delta t$ [fs]')
ylabel('$Q(\Delta t)$ [fC]')
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])

end
end

if  strcmp('auto',get(gca,'XLimMode')) && strcmp('auto',get(gca,'YLimMode')) == 1
    axis tight
    disp('setting axis to tight')
end
grid on
set(gca,'fontsize', 15)
set(gca,'Box','on')

if exist('POPOUT') == 1
%Orig: LEFT = 0.06; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;
LEFT = 0.15; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
%InSet = [0.105    0.12    0.018         0.02]; %Fontsize 14
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
end
if exist('POPOUT') == 0
LEFT = 0.06; BOTTOM = 0.13; RIGHT = 0.35; TOP = 0.05;    
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
%InSet = [0.105    0.12    0.018         0.02]; %Fontsize 14
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
end

%THESISGRAPH2 = 1
if exist('THESISGRAPH2') == 1
TGPH = figure(3);
TGPH.Position = [550 270 920 520];
hold on
plot(F_0,abs(Q)/10^-15)
% plot([0,5*10^10],[1.60217662*10^-4, 1.60217662*10^-4],'r--')    
xlabel('Field Strength [$Vm^{-1}$]')
ylabel('Photoinduced Charge [fC]')    
TGPH.Position = [550 270 920 540]
legend('$N_{cyc}$ = 1.0','$N_{cyc}$ = 1.5','$N_{cyc}$ = 2.0','$N_{cyc}$ = 2.5','$N_{cyc}$ = 3.0',...
       '$N_{cyc}$ = 3.5','$N_{cyc}$ = 4.0','$N_{cyc}$ = 4.5','$N_{cyc}$ = 5.0','$N_{cyc}$ = 5.5','$N_{cyc}$ = 6.0','e')
grid on
set(gca,'YScale','log')
set(gca,'fontsize', 15)
set(gca, 'FontName', 'SMU Serif')
set(gca,'Box','on')
LEFT = 0.12; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;    
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
end
%FIG4PLS = 1
if exist('FIG4PLS') == 1
TGPH = figure(4)
hold on

plt1 = plot(ETI.ttt.tdisp/10^-15,real(circshift(Ediresdt,ETD.tcirc)),'c');%,'FaceColor','b','LineStyle','none' %Total overlap = area
plt2 = plot(ETD.ttt.tdisp/10^-15,real(circshift(ETD.ttt.Etdshift,ETD.tcirc)),'r','LineWidth',1); %No Overlap - Shifted pulse
plt3 = plot(ETI.ttt.tdisp/10^-15,real(0.2*ETI.ttt.Etdispc),'g','LineWidth',1); %No Overlap - Stationary pulse
legend('Sum','Fundamental','$2_{nd}$ Harmonic')
%plt4 = plot(Delt./(N_Delt/N),QHarm/max(QHarm)+1,'m'); %Overlap
%plot(ETI.ttt.tdisp,real(Ediresdt),'b','LineWidth',1.3);

xlabel('$\Delta t$ [fs]')
ylabel('Sum of Vector Potential')
axis([ETI.ttt.tdisp(3400)/10^-15 ETI.ttt.tdisp(end-3400)/10^-15 -1.2 2])

%xlabel('\Delta t [s]')
TGPH.Position = [550 270 560 520]
grid on
set(gca,'fontsize', 16)
set(gca, 'FontName', 'Computer Modern')
set(gca,'Box','on')
LEFT = 0.15; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;    
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
end
