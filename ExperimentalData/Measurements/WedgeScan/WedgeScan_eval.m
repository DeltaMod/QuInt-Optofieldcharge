clear all;clc;
%
% add global scripts
addpath('C:\Users\fa5471la\Local_Matlab\global_scripts');
bd = '\\fysfile01\atomhome$\fa5471la\My Documents\data\2019_04_26_wedgescan\';
bd2 = 'C:\Users\fa5471la\Local_Matlab\2019_04_26_LWC-WedgeScan\';
%
%
%M=dlmread([bd 'test_lowamp.dat']);
M2=dlmread([bd 'test_lowamp_nolock.dat']);
%
%glass = M(:,1);
%signal = -M(:,2); %amps/1000 or so
glass2 = M2(:,1);
signal2 = -M2(:,2); %amps/1000 or so
%
figure(1);clf;
    plot(glass,signal,'-b'); hold on;
    plot(glass2,signal2,'-r');
%
figure(2);clf;
    plot(glass,signal-signal2);
% fourier
[nu,gl_nu] = FFT(glass,signal-mean(signal));
[nu,gl2_nu] = FFT(glass,signal2-mean(signal2));
figure(3);clf;
    plot(nu,abs(gl_nu)); hold on;
    plot(nu,abs(gl2_nu));
    plot([1/50 1/50],[-1e3 1e4],'--k')
    ylim([0 7000])
% bp
sig_filt = real(bp(glass,signal,[0.005 0.065],1e-3));
sig2_filt = real(bp(glass,signal2,[0.005 0.065],1e-3));
figure(4);clf;
    plot(glass,signal,'-b'); hold on;
    plot(glass,sig_filt,'-r');
    plot(glass,sig2_filt+2e4,'-k')
%%
M=dlmread([bd 'longscan.dat']);
M2=dlmread([bd 'longscan_nolock.dat']);
%
glass = M(:,1);
signal = -M(:,2); %amps/1000 or so
glass2 = M2(:,1);
signal2 = -M2(:,2); %amps/1000 or so
%
figure(1);clf;
    plot(glass,signal,'-b'); hold on;
    plot(glass2,signal2,'-r');
    xlabel('Glass insertion (µm)')
    ylabel('Signal (arb.u.)')
    legend('w. locking','w.o. locking')
    xlim([-1.1e3 1.1e3])
saveFigure([bd2 'longscan_w-wo-lock.png'],[14 10])
%
figure(2);clf;
    plot(glass,signal-signal2);
% fourier
[nu,gl_nu] = FFT(glass,signal-mean(signal));
[nu,gl2_nu] = FFT(glass,signal2-mean(signal2));
figure(3);clf;
    semilogy(nu,abs(gl_nu),'-b'); hold on;
    plot(nu,abs(gl2_nu),'-r');
    ylim([1e-10 1e-7])
    xlim([0 0.1])
%     plot([1/50 1/50],[-1e3 1e4],'--k')
%     ylim([0 7000])
% bp
sig_filt = real(bp(glass,signal,[0.005 0.065],1e-3));
sig2_filt = real(bp(glass,signal2,[0.005 0.065],1e-3));
figure(4);clf;
    plot(glass,signal,'-b'); hold on;
    plot(glass,sig_filt,'-r');
    plot(glass,sig2_filt+3e-7,'-k')
figure(5);clf;
    plot(glass,sig_filt-sig2_filt)
%%
M=dlmread([bd 'test6.dat']);
% M2=dlmread([bd 'longscan_nolock.dat']);
%
glass = M(:,1);
signal = -M(:,2); %amps/1000 or so
% glass2 = M2(:,1);
% signal2 = -M2(:,2); %amps/1000 or so
%
figure(1);clf;
    plot(glass,signal,'-b'); hold on;
%     plot(glass2,signal2,'-r');
%
% figure(2);clf;
%     plot(glass,signal-signal2);
% fourier
[nu,gl_nu] = FFT(glass,signal-mean(signal));
% [nu,gl2_nu] = FFT(glass,signal2-mean(signal2));
figure(3);clf;
    plot(nu,abs(gl_nu)); hold on;
%     plot(nu,abs(gl2_nu));
%     plot([1/50 1/50],[-1e3 1e4],'--k')
%     ylim([0 7000])
% bp
sig_filt = real(bp(glass,signal,[0.005 0.065],1e-3));
% sig2_filt = real(bp(glass,signal2,[0.005 0.065],1e-3));
figure(4);clf;
    plot(glass,signal,'-b'); hold on;
    plot(glass,sig_filt,'-r');
%     plot(glass,sig2_filt+3e-7,'-k')