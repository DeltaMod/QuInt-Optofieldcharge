%% This UI Document is used to plot the graphs, and nothing else! %%
if exist('POPOUT','var') == 1
FIG = figure(1); FIG.Renderer = 'painters'; FIG.Position = [480 270 960 540]; FIG.Name = PlotSelect{1};
set(gca,'fontsize', 12)
set(gca, 'FontName', 'Computer Modern')

box on
end
%% PhotoinducedChargetoUI %%
if SimSelect == "Simple Photoinduced Charge"
if exist('UIRUN','var') == 0; PlotSelect = "Gaussian Laser Pulse"; FIG = figure(1);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,1); end
if PlotSelect == "Gaussian Laser Pulse"
cla reset
plot((t-t_0)/T_x,a2n(:,1:N),'LineWidth',1.3)
ylim([-1 1])
set(gca, 'XLimSpec', 'Tight');
xlabel('Optical Cycle ((t-t_0)/T)')
ylabel('Vector Potential (and a^2n+1)')
legend('a(t)^3','a(t)^5','a(t)^7','a(t)^9','Location','best' )
grid on 
end

if exist('UIRUN','var') == 0; PlotSelect = "<a^2n+1>"; FIG = figure(2);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,3); end
if PlotSelect == "<a^2n+1>"
cla reset;
semilogy(Ncyc2,Vpota(:,:))
hold on
%axis([1 2.5 10^-6 10^0])
grid on
legend('<a^3>','<a^5>','<a^7>','<a^9>','Location','best' )
xlabel('Number of Optical Cycles (t/T)')
ylabel('Vector Potential momenta')
end

if exist('UIRUN','var') == 0; PlotSelect = "Photoinduced Charge";  FIG = figure(3); clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(2,3); end
if PlotSelect == "Photoinduced Charge"
    cla reset
plot(F_0,Q);
grid on
xlabel('Optical Field [Vm^-^1]')
ylabel('Photinduced Charge [C]')
end
end



%% DispersionToUI %%
if SimSelect == "Dispersion and Photoinduced Charge"
if exist('UIRUN','var') == 0; PlotSelect = "Non Fourier Et(t)"; FIG = figure(1);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(1,3); end
if PlotSelect == "Non Fourier Et(t)"
cla reset;
%subplot(1,2,1)
plot(t-t_0,real(Et)); xlabel('time [s]'); ylabel('E_t(t)'); fprintf(['\n','Plotting [',PlotSelect{1},']']);
axis tight
grid on
end

if exist('UIRUN','var') == 0; PlotSelect = "E_t FFT Plot"; FIG = figure(2);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(2,3); end
if PlotSelect == "E_t FFT Plot"
cla reset;
plot(Estr.ttt.w,abs(Estr.ttt.Ew)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('\omega [rad/s]')
ylabel('Amplitude')
hold on
plot([w_0x,w_0x],[0,abs(max(Estr.ttt.Ew))])
end

if exist('UIRUN','var') == 0; PlotSelect = "E_t Final Applied Dispersion"; FIG = figure(3);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end 
if PlotSelect == "E_t Final Applied Dispersion"
cla reset;
plot(linspace(t1,t2,length(Estr.ttt.tdisp))-t_0,real(Estr.ttt.Etdisp));grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('s[fs]')
ylabel('amplitude')
grid on
end

if exist('UIRUN','var') == 0; PlotSelect = "E_t(t)^2n+1 After Dispersion"; FIG = figure(4);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,1); end
if PlotSelect == "E_t(t)^2n+1 After Dispersion"
cla reset
plot(Estr.ttt.tdisp,Etf2n(:,:),'LineWidth',1.3)
hold on
ylim([-1 1])
xlabel('Optical Cycle ((t-t_0)/T)')
ylabel('Vector Potential (and a^2n+1)')
legend('E_t(t)','E_t(t)^3','E_t(t)^5','E_t(t)^7','E_t(t)^9','Location','best' )
grid on 
end

if exist('UIRUN','var') == 0; PlotSelect = "Post Dispersion Photoinduced Charge";  FIG = figure(5); clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(2,3); end
if PlotSelect == "Post Dispersion Photoinduced Charge"
    cla reset
plot(F_0,abs(Q));
grid on
xlabel('Optical Field [Vm^-^1]')
ylabel('Photinduced Charge [C]')
end

if exist('UIRUN','var') == 0; PlotSelect = "E_w IFFT Dispersion"; FIG = figure(6);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(1,1); end 
if PlotSelect ==  "E_w IFFT Dispersion"
cla reset;
plot(Estw.ftt.tdisp,real(Estw.ftt.Etdisp)); grid on; fprintf(['\n','Plotting [',PlotSelect{1},']']);
legend(['\phi_2 = ' num2str(SiO2.phi(3),'%.4g')])
xlabel('Time [s]')
ylabel('Amplitude')
grid on
end


if exist('UIRUN','var') == 0; PlotSelect = "Current Pulse Overlap at given Delta t"; FIG = figure(7);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end
if PlotSelect == "Current Pulse Overlap at given Delta t"
SldInd = find(ismember(Deltn,Slidern));
ETD = E_td; ETI = E_ti;
ETD.ttt.Etdshift = circshift(ETD.ttt.Etdisp,Deltn(SldInd)); ETD.ttt.tdshift = circshift(ETD.ttt.tdisp,Deltn(SldInd)); 
Ediresdt = ETD.ttt.Etdshift+ETI.ttt.Etdisp;

for m = 1:4
    a2harmdt(m) = w_0x*trapz(ETD.ttt.tdisp,real(Ediresdt).^(2*m+1));
    a2poldt(m)  = w_0x*trapz(ETD.ttt.tdisp,real(ETD.ttt.Etdshift).*real(ETI.ttt.Etdisp).^(2*m)); 
end
QHarmdt(:) = fun_Q(F_0(end),F_a,a2harmdt(:),Aeff);
QPoldt(:)  = fun_QDelt(F_0(end)/12.5,F_0(end)/1.25,a2poldt(:),Aeff);
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
xlabel('\Delta t [s]')
ylabel('Sum of Vector Potential')
grid on
axis([ETI.ttt.tdisp(1) ETI.ttt.tdisp(end) -3 3])
hold off
end


if exist('UIRUN','var') == 0; PlotSelect = "Temporal Overlap Induced Charge"; FIG = figure(8);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end
if PlotSelect == "Temporal Overlap Induced Charge"
cla reset;
plot(Delt,real(QHarm));fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('\Delta t [s]')
ylabel('Q(\Delta t [C])')
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])
grid on
axis tight
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
xlabel('\Delta t [s]')
ylabel('Vector Potential Momenta')
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])
grid on
axis tight
%figure(1)
%plot(1,1)
end

if exist('UIRUN','var') == 0; PlotSelect = "Dual Polarised Induced Charge"; FIG = figure(10);clf; FIG.Name = PlotSelect{1}; FIG.Position = FigPos(3,2); end
if PlotSelect == "Dual Polarised Induced Charge"
cla reset;
plot(Delt,real(QPol));fprintf(['\n','Plotting [',PlotSelect{1},']']);
xlabel('\Delta t [s]')
ylabel('Q(\Delta t) [C]')
%legend(['T_E_i =' num2str((0.5*n-10)*T+T,'%.4g') ')'])
grid on
axis tight
end
end