POW  =   [0 6 7 8 9 10 11 12 13 14 15 16 17 18 19.2] 
VOLT = [0 0.27 0.48 0.74 1.37 1.9 3.1 4.6 6.7 8.8 12.6 15.9 20.0 23.4 26.6]

FIG = figure(1); FIG.Renderer = 'painters'; FIG.Position = [480 270 960 540]; clf ;
set(gca,'fontsize', 12)
set(gca, 'FontName', 'Computer Modern')

plot(POW,VOLT)
xlabel('Power [mW]')
ylabel('Voltage [mV]')
grid on
box on
hold on
axis tight

% %% Fitting Curves
% [Fit,Err] = polyfit([POW],[VOLT],3);
% X = linspace(0,25,1000)
% [Val] = polyval(Fit,X,Err);
% plot(X,Val,'r--')
% %Calculating the standard Deviation
% b_err = sqrt(diag((Err.R)\inv(Err.R'))./Err.normr.^2./Err.df);