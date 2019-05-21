close all;
addpath(genpath(fileparts(which('IVCurvesandResistance.m'))))
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex'); set(0, 'defaultTextInterpreter','latex');
set(0,'defaultAxesFontName', 'CMU Serif'); set(0,'defaultTextFontName', 'CMU Serif');

d    = dir('Measurements/DeviceIVs/*.txt');  % get the list of files
AMP = [1e-3, 1e-3, 1e-3, 500e-6, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3,200e-6, 200e-6, 1e-3, 500e-6, 500e-6]
for n = 1:length(d)
FN{n} = d(n).name;
if FN{n}(2) == '1'
    DN{n} = [FN{n}(1:3),FN{n}(end-7:end-4)]
    
else
DN{n} = [FN{n}(1:2),FN{n}(end-7:end-4)];
end
Dat.(DN{n}).Amp = AMP(n);
end



   

for n = 1:length(DN)
fid = fopen(FN{n},'rt');
tLines = fgetl(fid); tLines = fgetl(fid);
numCols = numel(strfind(tLines,sprintf('\t'))) + 1; fclose(fid);

fid = fopen(FN{n},'rt');
indata = textscan(fid, '%f', 'HeaderLines',2,'Delimiter','\t');

fclose(fid);    

Csplt = linspace(1,length(indata{1})-numCols+1,(length(indata{1}))/numCols);
        Dat.(DN{n}).t = indata{1}(Csplt); 
        Dat.(DN{n}).V = indata{1}(Csplt+1);
        Dat.(DN{n}).I1 = indata{1}(Csplt+2);
        Dat.(DN{n}).I2 = indata{1}(Csplt+3);
end




figure(1)
hold on
for n = 1:length(DN)
L = length(Dat.(DN{n}).V);
R = round(L/4):1:round(L-L/4);
[Fit,Err] = polyfit(Dat.(DN{n}).V(R),Dat.(DN{n}).I1(R),1);
X = linspace(Dat.(DN{n}).V(1),Dat.(DN{n}).V(end),50);
[Val] = polyval(Fit,X,Err);
plot(X,Val,'r')
%Calculating the standard Deviation
b_err = sqrt(diag((Err.R)\inv(Err.R'))./Err.normr.^2./Err.df);
Dat.(DN{n}).R = 1/Fit(1);
Dat.(DN{n}).Err = b_err;
end
hold on
grid on
for n = 1:length(DN)
    plot(Dat.(DN{n}).V,Dat.(DN{n}).I1)
end




figure(3)
hold on
grid on
for n = 1:length(DN)
    %bar(n,Dat.(DN{n}).R)
    bar(n,Dat.(DN{n}).R)
end
    TABLE = ['Device',' & ', 'Resistance [\Omega] \\'];

for n = 1:length(DN)
    TABLE = [TABLE,DN{n},' & ',num2str(round(Dat.(DN{n}).R,5)),' \\'];
end
disp(TABLE)
  

figure(4)

plot(Dat.(DN{2}).V,Dat.(DN{2}).I1*10^3)
hold on
%plot(Dat.(DN{3}).V,Dat.(DN{3}).I1)
plot(Dat.(DN{4}).V,Dat.(DN{4}).I1*10^3)
plot(Dat.(DN{5}).V,Dat.(DN{5}).I1*10^3)
plot(Dat.(DN{6}).V,Dat.(DN{6}).I1*10^3)
LEG = legend([DN(2),DN(4:6)],'Location','best')
title(LEG,'Device Names')
xlabel('Voltage [V]'); ylabel('Current [mA]');
set(gca,'fontsize', 16)
set(gca, 'FontName', 'Computer Modern')
set(gca,'Box','on')
grid on
axis tight

LEFT = 0.12; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;    
InSet = [LEFT BOTTOM RIGHT TOP]; %Fontsize 12
%InSet = [0.105    0.12    0.018         0.02]; %Fontsize 14
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])