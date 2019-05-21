function Ef = bp(t, E, bp_nu, cutWidth)
% >> Ef = bp(t, E, bp_nu, cutWidth)
% Bandpass filter
%
%   t = a 1xN vector (e.g. time in s)
%   E = a 1xN vector (electric field of transient)
%   bp_nu = a 1x2 vector with the cut-on and cut-off frequency
%   cutWidth = width of exponential to smooth edges (default: 1)
%   
% Return values:
%
%   Ef = 1xN vector with the filtered data

% Author: Olaf Schubert (e-Mail: olaf.schubert@uni-konstanz.de)
% Date: 15.07.2010
%
% Version history (add to this list if you fix bugs or change something!):
% -----------------------------
% 1.0.0    first version
%

%
global ENABLEDEBUG;

if nargin<4
    cutWidth = 1;
end
%% Interpolate data and do a fourier transformation
    [nu En ti] = fftO(t, E, length(E));
% create filter mask
    mask = createMask( nu, bp_nu, cutWidth);
% filter data
    En2 = En.*mask;
% interpolate filtered data on old scale
    Ef = interp1(ti,ifft(ifftshift(En2)),t);
    
    if ENABLEDEBUG>0
       figure(99);clf;
        plot(nu,mask);
    end