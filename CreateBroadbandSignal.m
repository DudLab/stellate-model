function [f_y,t] = CreateBroadbandSignal(f_max,duration,dt);
% 
% This is the most simple version of a broadband signal generator where you
% take random white noise and filter it down to broadband signals. Future
% development will aim to create more interesting filters that can be
% used, but for the moment simple brickwall filters are used.
% 
% Time variable are expressed in ms (for consistency) and f_max is
% expressed in Hz

% Determine total points
t           = dt:dt:duration;
totalp      = length(t);

% Set the FFT length (2048 gives a max f of 2 kHz)
fft_n       = totalp;

% Get sampling rate
fixrate     = 1000/dt;

s_y         = randn(1,totalp);
S_y         = fft(s_y,fft_n);
F_y         = zeros(1,totalp);
F_y(1:round(((f_max/dt))))  = S_y(1:round(((f_max/dt))));
f_y         = real(ifft(F_y,fft_n));
f_y         = f_y ./ max(abs(f_y),[],2);

(fixrate/fft_n)*f_max

figure(1);
plot(t,f_y);
drawnow;