close all; clear all;

addpath(genpath(fileparts(which('IVCurvesandResistance.m'))))
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex'); set(0, 'defaultTextInterpreter','latex');
set(0,'defaultAxesFontName', 'CMU Serif'); set(0,'defaultTextFontName', 'CMU Serif');

d    = dir('Measurements/Omega2Omega_DelayScan/*.txt');  % get the list of files

for n = 1:length(d)
FN{n} = d(n).name;
DN{n} = [FN{n}(1:7),num2str(n)];
end


for n = 1:length(DN)
fid = fopen(FN{n},'rt');
tLines = fgetl(fid); tLines = fgetl(fid);
numCols = numel(strfind(tLines,sprintf('\t'))) + 1; fclose(fid);

fid = fopen(FN{n},'rt');
indata = textscan(fid, '%f', 'HeaderLines',0,'Delimiter','\t');

fclose(fid);    

Csplt = linspace(1,length(indata{1})-numCols+1,(length(indata{1}))/numCols);
        Dat.(DN{n}).V = indata{1}(Csplt); 
        Dat.(DN{n}).I = indata{1}(Csplt+1);
end

figure(1)
hold on
for n = 2
plot(Dat.(DN{n}).V,Dat.(DN{n}).I/-10^-9)
end
grid on


LEFT = 0.1; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;
InSet = [LEFT BOTTOM RIGHT TOP];
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
grid on
set(gca,'fontsize', 16)
set(gca, 'FontName', 'Computer Modern')
set(gca,'Box','on')
axis tight
xlabel('Mirror Position [mm]')
ylabel('Current [nA]')

clear all

d    = dir('Measurements/WedgeScan/*.dat');  % get the list of files
for n = 1:length(d)
FN{n} = d(n).name;
DN{n} = FN{n}(1:end-4);
end


for n = 1:length(DN)
fid = fopen(FN{n},'rt');
tLines = fgetl(fid); tLines = fgetl(fid);
numCols = numel(strfind(tLines,sprintf('\t'))) + 1; fclose(fid);

fid = fopen(FN{n},'rt');
indata = textscan(fid, '%f', 'HeaderLines',0,'Delimiter','\t');

fclose(fid);    

Csplt = linspace(1,length(indata{1})-numCols+1,(length(indata{1}))/numCols);
        Dat.(DN{n}).V = indata{1}(Csplt); 
        Dat.(DN{n}).I = indata{1}(Csplt+1);
end

figure(2)
hold on
for n = 3:4
    hold on
plot(Dat.(DN{n}).V,Dat.(DN{n}).I/-10^-9)
end
grid on
LEFT = 0.1; BOTTOM = 0.13; RIGHT = 0.05; TOP = 0.05;
InSet = [LEFT BOTTOM RIGHT TOP];
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
grid on
set(gca,'fontsize', 16)
set(gca, 'FontName', 'Computer Modern')
set(gca,'Box','on')
axis tight
legend('With Lock-in', 'Without Lock-in','Location','Best')
xlabel('Glass Insertion [$\mu m$]')
ylabel('Current [nA]')