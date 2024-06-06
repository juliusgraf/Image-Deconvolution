%-------------------------------------------------------------------------
%----- Fichier : deconv_naive
%----- Objet   : Deconvolution naive
%----- References principales : Idier, 
%-------------------------------------------------------------------------

close all;
clear all;

L1x = 100;	    % duree d'observation du signal d'entree (en s)
fe = 1;		    % frequence d'echantillonnage (Hz)
Te = 1/fe;		% periode d'echantillonnage
N1x = L1x/Te;	% nombre de points du signal x
kx = 0:N1x;	    % vecteur des indices temporels  
tx= kx*Te;		% vecteur des instants d'echantillonnage


% 1. Reponse impulsionnelle h d'un filtre gaussien
%
mu_h=15*Te;
% la largeur de la reponse impulsionnelle va rendre le probleme plus ou
% moins difficile
sigma_h=5*Te;
tt=(0:Te:30*Te);
h=(1/(sigma_h*sqrt(2*pi)))*exp(-(((tt-mu_h)/(sqrt(2)*sigma_h)).*((tt-mu_h)/(sqrt(2)*sigma_h))));

figure(1);
subplot(2,3,2)
plot(tt,h)
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Reponse impulsionnelle $h(t)$', 'Interpreter', 'latex');

% 2. Signal x
%
k1 = 0:N1x/2-1;

x1 = zeros(1,length(k1));
k2 = N1x/2:N1x;  % modif  N1x-1 => N1x
x2 = ones(1,length(k2));

x = [x1 x2];
subplot(2,3,1)
plot(tx,x);
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal entree $x(t)$', 'Interpreter', 'latex');

%3. Convolution
%

y_nb=conv(x,h);
N1=length(y_nb);
k=0:1:N1-1;
t=k*Te;
subplot(2,3,3)
plot(t,y_nb)
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal sortie non bruite $y_{nb}(t)$', 'Interpreter', 'latex');

% 4. Bruitage
%
SNR = 30;
y = adgnoise(y_nb, SNR);
subplot(2,3,6)
plot(t,y);
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal sortie bruite $y(t)$', 'Interpreter', 'latex');

% 5. Representation du bruit 
%
w = y - y_nb;
subplot(2,3,5)
plot(t,w); 
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Bruit $w(t)$', 'Interpreter', 'latex');
RSB=10*log10((y_nb*y_nb'/length(y_nb))/var(w)) % 



% 6. Reponse en frequence du filtre
%
N=1024;
n = -N/2:N/2-1;
f = n*fe/N;
H = fft(h,N);
figure(2); subplot(2,3,2)
plot(f,fftshift(20*log10(abs(H))));
xlabel('frequence (Hz)', 'Interpreter', 'latex'); ylabel('dB', 'Interpreter', 'latex'); title('dse du filtre $H(f)$', 'Interpreter', 'latex');


% 7. TFD X du signal d'entree, TFD Y du signal de sortie, TFD X_rec du signal reconstruit 
%
X = fft(x,N);
subplot(2,3,1)
plot(f,fftshift(10*log10(abs(X).^2/length(x))));
xlabel('frequence (Hz)', 'Interpreter', 'latex'); ylabel('dB', 'Interpreter', 'latex'); title('dsp $\Gamma_x(f)$', 'Interpreter', 'latex');

Y = fft(y,N); % spectre du signal bruite
%%%Y=fft(y_nb,N); % spectre du signal non bruite
subplot(2,3,3)
plot(f,fftshift(10*log10(abs(Y).^2/length(y))));
xlabel('frequence (Hz)', 'Interpreter', 'latex'); ylabel('dB', 'Interpreter', 'latex'); title('dsp $\Gamma_y(f)$', 'Interpreter', 'latex');


% Reponse en frequence du filtre inverse
%
Hinv = 1./H;
hinv = real(ifft(Hinv));
subplot(2,3,5)
plot(f,fftshift(20*log10(abs(Hinv))));
xlabel('frequence (Hz)', 'Interpreter', 'latex'); ylabel('dB', 'Interpreter', 'latex'); title('dse du filtre inverse $1/H$', 'Interpreter', 'latex');

% DSP du bruit
%
W = fft(w,N);
subplot(2,3,6)
plot(f,fftshift(10*log10(abs(W).^2/length(w))));
xlabel('frequence (Hz)', 'Interpreter', 'latex'); ylabel('dB', 'Interpreter', 'latex'); title('dsp $\Gamma_w(f)$', 'Interpreter', 'latex');

% calcul du signal reconstruit par filtrage inverse
%
X_rec2 = Y./H;
x_rec2 = real(ifft(X_rec2));

Nrec2=length(x_rec2);
krec2=0:1:Nrec2-1;
trec2=krec2*Te;
figure(1);
subplot(2,3,4)
plot(trec2(1:length(x)),x_rec2(1:length(x)))
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal d''entree reconstruit par filtrage inverse', 'Interpreter', 'latex');



figure(2)
% subplot(2,3,3)
% plot(f,fftshift(10*log10(abs(Y).^2/length(y))));
% xlabel('frequence (Hz)'); ylabel('dB'); title('dsp Y(f)');
% 
% 
% subplot(2,3,1)
% plot(f,fftshift(10*log10(abs(X.^2)/length(x))));
% xlabel('frequence (Hz)'); ylabel('dB'); title('dsp X(f)');
% 
% subplot(2,3,2)
% plot(f,fftshift(20*log10(abs(H))));
% xlabel('frequence (Hz)'); ylabel('dB'); title('dse du filtre H(f)');
% 
% subplot(2,3,6)
% plot(f,fftshift(10*log10(abs(W).^2/length(w))));
% xlabel('frequence (Hz)'); ylabel('dB'); title('dsp W(f)');
% 
% subplot(2,3,5)
% plot(f,fftshift(20*log10(abs(Hinv))));
% xlabel('frequence (Hz)'); ylabel('dB'); title('dse du filtre inverse 1./H');
% 
subplot(2,3,4)
plot(f,fftshift(10*log10(abs(X_rec2).^2/length(x_rec2))));
xlabel('frequence (Hz)', 'Interpreter', 'latex'); ylabel('dB', 'Interpreter', 'latex'); title('dsp $\Gamma_{x_{rec2}}(f)$', 'Interpreter', 'latex');
