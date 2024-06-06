close all;
clear all;

L1x = 100;	    % duration of observation of input signal (in sec.)
fe = 1;		    % sampling frequency (Hz)
Te = 1/fe;		% sampling period
N1x = L1x/Te;	% number of points in signal x
kx = 0:N1x;	    % time indices  
tx= kx*Te;		% sampling times 


% 1. Impulse response h of a Gaussian filter
%
mu_h=15*Te;
% the width of impulse response makes the problem more or less difficult
sigma_h=5*Te;
tt=(0:Te:30*Te);
h=(1/(sigma_h*sqrt(2*pi)))*exp(-(((tt-mu_h)/(sqrt(2)*sigma_h)).*((tt-mu_h)/(sqrt(2)*sigma_h))));

figure(1);
subplot(2,3,2)
plot(tt,h)
xlabel('time (s)'); ylabel('amplitude'); title('Impulse Response h(t)');

%% 2. Signal x
%

k1 = 0:N1x/2-1;
x1 = zeros(1,length(k1));
k2 = N1x/2:N1x;  % modif  N1x-1 => N1x
x2 = ones(1,length(k2));

x = [x1 x2];
subplot(2,3,1)
plot(tx,x);
xlabel('time (s)'); ylabel('amplitude'); title('Imput signal x(t)');

%% 3. Convolution
%

y_nb=conv(x,h,'full');
y_nb=y_nb(end-100:end);
N1=length(y_nb);
k=0:1:N1-1;
t=k*Te;
subplot(2,3,3)
plot(t,y_nb)
xlabel('time (s)'); ylabel('amplitude'); title('Noiseless output signal y_{nb}(t)');

%% 4. Noise addition
%
SNR = 30;
[y,sigma_br] = adgnoise(y_nb,SNR);
  
subplot(2,3,6)
plot(t,y);
xlabel('time (s)'); ylabel('amplitude'); title('Noisy signal y(t)');

%% 5. Frequency response of filter
%
N= 101;
f = -N/2:N/2-1;
H = fft(h,N);
figure(2); subplot(2,3,2)
plot(f,fftshift(20*log10(abs(H))));
xlabel('frequency (Hz)'); ylabel('dB'); title('ESD of filter H(f)');

%% 6. Frequency response of data
%
Y = fft(y,N);
subplot(2,3,1)
plot(f,fftshift(20*log10(abs(H))));
xlabel('frequency (Hz)'); ylabel('dB'); title('ESD of filter Y(f)');

%% Frequency response of reconstructed signal
X_rec = Y./H;
x_rec = real(ifft(X_rec));
Nrec=length(x_rec);
krec=0:1:Nrec-1;
trec=krec*Te;
figure;
plot(x_rec(1:length(x))); 
xlabel('time (s)'); ylabel('amplitude'); title('Reconstruced input signal by inverse filtering');

%% Part 2 bis

d = [1;-1];
D = fft(d,101)';
alpha = 10^(-2.2);
g = H ./ (abs(H).^2 + alpha * abs(D).^2);
X_rec2 = Y.*g;
x_rec2 = real(ifft(X_rec2));

Nrec2=length(x_rec2);
krec2=0:1:Nrec2-1;
trec2=krec2*Te;
figure;
plot(trec2(1:length(x)),x_rec2(1:length(x)))
xlabel('time (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Reconstruced input signal by inverse filtering with $\alpha=1$', 'Interpreter', 'latex');
