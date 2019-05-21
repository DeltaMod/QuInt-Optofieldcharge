%% This function performs a FFT of a e.g. time depending signal into the frequency space. The sampling rate is given
%  over the data interval 

function [xout,funcout]=FFT(xvalues,yvalues)
M=length(xvalues);
funcout = fftshift((fft(yvalues))./M);
Fsampling=1/((xvalues(M)-xvalues(1))./M);
xout = .5*abs(Fsampling).*linspace(-1,1,M);

% % nextpow2 code
% M=length(xvalues);
% NFFT = 2.^nextpow2(M);
% funcout = fftshift((fft(yvalues,NFFT)/M));
% Fsampling=1/((xvalues(M)-xvalues(1))./M);
% xout = .5*abs(Fsampling).*linspace(-1,1,NFFT);
end
