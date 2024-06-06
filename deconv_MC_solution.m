%-------------------------------------------------------------------------
%----- Fichier : deconv_MC
%----- Objet   : Deconvolution par moindres carres regularises ou non
%----- References principales : Idier
%-------------------------------------------------------------------------

close all;
clear all;

L_x = 100;	    % duree d'observation du signal d'entree
fe = 1;		    % frequence d'echantillonage
Te = 1/fe;		% periode d'echantillonage
N_x = L_x/Te;	% nombre de points du signal x
k_x = 0:N_x;	% index temporel
t_x= k_x*Te;    % base de temps
%**************************************************************************
% 1. Signal x
k1 = 0:N_x/2-1;
x1 = zeros(1,N_x/2);
k2 = N_x/2:N_x;
x2 = ones(1,length(k2));
x = [x1'; x2'];
subplot(2,3,1)
plot(t_x,x);
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal entree $x(t)$', 'Interpreter', 'latex');


%**************************************************************************
% 2. Reponse impulsionnelle h de filtres lineaires gaussiens
% test de plusieurs filtres gaussien ; la largeur de la reponse impulsionnelle va rendre le probleme plus ou
% moins difficile
%**************************************************************************
N_sigma=30;
mu_h=15*Te;
for sigma=1:N_sigma,
   sigma_h=sigma*Te;
   L_h=30*Te;
   t_h=(0:Te:L_h);
   N_h=length(t_h);
   h=(1/(sigma_h*sqrt(2*pi)))*exp(-(((t_h-mu_h)/(sqrt(2)*sigma_h)).*((t_h-mu_h)/(sqrt(2)*sigma_h))));
   
   %-----------------------------------------------------------------------
   % 3. Construction de la matrice H de convolution (Toeplitz)
   N_x=length(x);
   N_y = N_x + N_h -1;
   
   hcol_1 = zeros(1,N_y);
   hcol_1(1:N_h) = hcol_1(1:N_h) + h;
   hlig_1 = zeros(1,N_x);  
   hlig_1(1,1) = hcol_1(1,1);
   H = toeplitz(hcol_1,hlig_1);
   cond_H(sigma)=cond(H);
end

%**************************************************************************
% Construction d'un seul filtre gaussien
%**************************************************************************
sigma=4;
sigma_h=sigma*Te;
L_h=30*Te;
t_h=(0:Te:L_h);
N_h=length(t_h);
h=(1/(sigma_h*sqrt(2*pi)))*exp(-(((t_h-mu_h)/(sqrt(2)*sigma_h)).*((t_h-mu_h)/(sqrt(2)*sigma_h))));
subplot(2,3,2)
plot(t_h,h)
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Reponse impulsionnelle $h(t)$', 'Interpreter', 'latex');

%--------------------------------------------------------------------------
%-
% 3. Construction de la matrice H de convolution (Toeplitz)

N_y = N_x + N_h -1;


hcol_1 = zeros(1,N_y);
hcol_1(1:N_h) = hcol_1(1:N_h) + h;
hlig_1 = zeros(1,N_x);  
hlig_1(1,1) = hcol_1(1,1);
H = toeplitz(hcol_1,hlig_1);
cond_H(sigma)=cond(H);

%**************************************************************************
%4. Convolution sous forme matricielle y_nb = H x

y_nb=H*x;
k=0:1:N_y-1;
t=k*Te;
subplot(2,3,3)
plot(t,y_nb)
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal sortie non bruite $y_{nb}(t)$', 'Interpreter', 'latex');

%**************************************************************************
% 5. Bruit additif gaussien

y= adgnoise(y_nb, 20);   %%% bruitage gaussien à 20 dB

% representation temporelle du signal de sortie bruite
subplot(2,3,6)
plot(t,y);
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal sortie bruite $y(t)$', 'Interpreter', 'latex');

%**************************************************************************
% 6. Representation du bruit 
w = y - y_nb;
RSB_y=10*log10((y'*y/N_y)/var(w));
subplot(2,3,5)
plot(t,w); 
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Bruit $w(t)$', 'Interpreter', 'latex');

%**************************************************************************
% 7. Deconvolution l2 par moindres carres

%--------------------------------------------------------------------------
x_rec=(H'*H)\(H'*y);
%--------------------------------------------------------------------------

subplot(2,3,4)
plot(t_x,x_rec)
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal d''entree reconstruit par MC', 'Interpreter', 'latex');

erreur_reconstruction = x - x_rec;
%print(erreur_reconstruction' * erreur_reconstruction);

%%
%**************************************************************************
% 8. Deconvolution par moindres carres regularises l2, penalisation des
% differences premieres

% Construction de la matrice D1 de differentiation

dcol_1 = zeros(1,N_x-1);
dcol_1(1) = 1; 
dlig_1 = zeros(1,N_x);  
dlig_1(1,1:2) = [1 -1];
D1 = toeplitz(dcol_1,dlig_1);


% Intervalle de variation du coefficient de regularisation
min_alpha=-7;
pas_alpha=0.1;
max_alpha=+2;
i_alpha=0;



% ligne suivante (ainsi que "end") à decommenter pour faire varier alpha

  for var_alpha=min_alpha:pas_alpha:max_alpha,

alpha=10^var_alpha;
i_alpha=i_alpha+1;
%--------------------------------------------------------------------------
x_rec_l2 (:,i_alpha)= (H'*H + alpha * D1'*D1)\(H'*y);
%--------------------------------------------------------------------------
err_rec(:,i_alpha)=x-x_rec_l2(:,i_alpha); % erreur de reconstruction

Werr_rec(i_alpha)=err_rec(:,i_alpha)'*err_rec(:,i_alpha); % energie de l'erreur de reconstruction
F(i_alpha)=norm(D1*x_rec_l2 (:,i_alpha),2)^2;
G(i_alpha)=norm(y-H*x_rec_l2 (:,i_alpha),2)^2;
end

% trace du signal reconstruit superpose au signal recherche
figure
clf,
  plot(t_x,x)
hold on
plot(t_x,x_rec_l2,'r')
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal d''entree reconstruit MC regularises', 'Interpreter', 'latex');
hold off
%pause;

%**************************************************************************
% energie de l'erreur de reconstruction en fonction du coefficient de regularisation
% à decommenter lorsque l'on fait varier alpha
figure
var_alpha=min_alpha:pas_alpha:max_alpha;
plot(var_alpha,10*log10(Werr_rec))
xlabel('log10$(\alpha)$', 'Interpreter', 'latex'); ylabel('energie erreur de reconstruction (dB)', 'Interpreter', 'latex'); title('Erreur de reconstruction en fonction de $\alpha$', 'Interpreter', 'latex');


%**************************************************************************
% trace du signal reconstruit superpose au signal recherche pour alpha
% optimal (choix "optimal" pour alpha = )
[W_err_rec_min,i_alpha_opt] = min(Werr_rec);
alpha_opt=10^(min_alpha+pas_alpha*(i_alpha_opt-1))
figure
plot(t_x,x)
hold on
plot(t_x,x_rec_l2(:,i_alpha_opt),'r')
xlabel('temps (s)', 'Interpreter', 'latex'); ylabel('amplitude', 'Interpreter', 'latex'); title('Signal d''entree reconstruit regularise pour $\alpha$ optimal', 'Interpreter', 'latex');
hold off

%**************************************************************************
% comparaison des conditionnements
figure
plot(log10(svd(H'*H)))
hold on
plot(log10(svd(H'*H + alpha_opt * D1'*D1)),'g')
xlabel('index'); ylabel('log10(valeurs singulieres)', 'Interpreter', 'latex'); 
title('valeurs singulieres', 'Interpreter', 'latex')
legend('$H^\top H$','$H^\top H + \alpha D_1^\top D_1$', 'Interpreter', 'latex');
hold off

display('conditionnement')
cond(H'*H)
%max(svd(H'*H))/min(svd(H'*H))
cond(H'*H + alpha_opt * D1'*D1)

%**************************************************************************
% Courbe en L
figure
plot(log10(G),log10(F),'.')
xlabel('log10(critere moindres carres)', 'Interpreter', 'latex'); ylabel('log10(terme de penalite)', 'Interpreter', 'latex'); title('Courbe en L', 'Interpreter', 'latex');



%**************************************************************************
% conditionnement en fonction de la largeur du filtre

figure
sigma=1:N_sigma;
stem(sigma,log10(cond_H))
grid
xlabel('$\sigma_h$', 'Interpreter', 'latex'); ylabel('log10(conditionnement de $H$)', 'Interpreter', 'latex'); title('Conditionnement en fonction de la largeur de la reponse impulsionnelle');

%%
%**************************************************************************
% 9. Approximation circulante et algorithmes rapides

N_x = length(x); % length of the input signal
Nh = length(h); % length of the impulse response
H = toeplitz([h zeros(1, N_x - Nh)], [h(1) zeros(1, N_x - 1)]);

% Fast inversion using FFT
lambda_h = fft(h, N_x)';
Y = fft(y, N_x);
X_hat = Y ./ lambda_h;
x_hat = real(ifft(X_hat));

% Plot the results
figure;
subplot(2,1,1);
plot(t_x, x);
title("Signal d'entree", 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_x, x_hat);
title('Signal reconstruit par FFT inverse', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Amplitude');

%**************************************************************************
% Moindres carrés régularisés
lambda_h = fft(h, N_x)';
% FFT de la dérivée première [1 -1]
d = [1 -1];
lambda_d = fft(d, N_x)';
% Calcul de gMCR
alpha = 0.1; % Exemple de valeur de régularisation
gMCR = lambda_h ./ (abs(lambda_h).^2 + alpha * abs(lambda_d).^2);
% Produit terme à terme des vecteurs Y et gMCR
X_hat = Y .* gMCR;
% Inverse FFT pour obtenir le signal reconstruit
x_hat = real(ifft(X_hat));

% Tracé des résultats
figure;
subplot(2,1,1);
plot(t_x, x);
title("Signal d'entree", 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_x, x_hat);
title('Signal reconstruit par regularisation FFT inverse pour $ \alpha =30$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Amplitude');

%%
%**************************************************************************
% 10. Déconvolution d’images par régularisation quadratique 

% Charger l'image cameraman.tif
x = imread('cameraman.tif');
x = im2double(x); % Convertir l'image en double précision

% Taille de l'image
[N, M] = size(x);

% Définir le filtre moyenneur 5x5
h1 = ones(5, 5) / 25;

% Filtrer l'image avec le filtre moyenneur
y_flou = conv2(x, h1, 'same');

% Ajouter du bruit additif blanc gaussien pour obtenir un RSB de 30 dB
rsb = 30; % en dB
Puissance_signal = var(y_flou(:)); % Puissance du signal
Puissance_bruit = Puissance_signal / (10^(rsb / 10)); % Puissance du bruit
y_bruite = adgnoise(y_flou, rsb); % Génération du bruit

% Visualiser les images
figure;
subplot(1, 3, 1);
imshow(x, []);
title('Image nette $\mathbf{x}$', 'Interpreter', 'latex');

subplot(1, 3, 2);
imshow(y_flou, []);
title('Image floue', 'Interpreter', 'latex');

subplot(1, 3, 3);
imshow(y_bruite, []);
title('Image floue et bruitee $\mathbf{y}$', 'Interpreter', 'latex');

% Calculer et visualiser le transfert en fréquence du filtre
H1 = fft2(h1, N, M);
H1_shifted = fftshift(H1); % Centrer le spectre

% Magnitude of the spectrum, adding a small constant to avoid log of zero
H1_magnitude = 20 * log10(abs(H1_shifted) + eps);

% Generate normalized frequencies
frequencies = (-N/2:N/2-1) / N;

% Create a frequency grid for x and y axes
[X, Y] = meshgrid(frequencies, frequencies);

% Plot the frequency response magnitude with normalized and centered axes
figure;
pcolor(X, Y, H1_magnitude);
shading interp;
colormap(jet);
colorbar;
xlabel('Normalized Spatial Frequency \nu_x', 'Interpreter', 'latex');
ylabel('Normalized Spatial Frequency \nu_y', 'Interpreter', 'latex');
title('Magnitude Spectrum of the Averaging Filter |H_1(\nu_x, \nu_y)|', 'Interpreter', 'latex');

% Another way to visualize the frequency response
figure;
imagesc(H1_magnitude);
colormap('jet');
colorbar;
title('Frequency Response of the Averaging Filter $20\log_{10}|H_1(\nu_x, \nu_y)|$', 'Interpreter', 'latex');
xlabel('Spatial Frequency $\nu_x$', 'Interpreter', 'latex');
ylabel('Spatial Frequency $\nu_y$', 'Interpreter', 'latex');

%%
%**************************************************************************
% 11. Simulation numérique du débruitage

% Dimensions de l'image
[N, M] = size(y_bruite);
% Fast inversion using FFT
lambda_h = fft2(h1, N, M);
Y = fft2(y_bruite, N, M);
X_hat = Y ./ lambda_h;
x_hat = real(ifft2(X_hat));

% Affichage du résultat de la déconvolution inverse
figure;
imagesc(x_hat)
colormap('gray')
title('Image reconstruite par déconvolution inverse');
set(gcf, 'Position', [100, 100, 800, 600]); % Redimensionner la fenêtre de la figure

% Définir le filtre de régularisation d1 (laplacien)
d1 = [0 -1 0; -1 4 -1; 0 -1 0];
lambda_d = fft2(d1, N, M); % FFT 2D de d1

% Paramètre de régularisation
alpha = 10;

% Calcul de gMCR en tenant compte de la régularisation
gMCR = lambda_h ./ (abs(lambda_h).^2 + alpha * abs(lambda_d).^2);

% Produit terme à terme des vecteurs Y et gMCR
X_hat_regularise = Y .* gMCR;

% Transformée inverse pour obtenir l'image déconvoluée
x_regularise = real(ifft2(X_hat_regularise));

% Affichage du résultat de la déconvolution régularisée
figure;
imagesc(x_regularise)
colormap('gray')
title(['Image reconstruite par régularisation FFT inverse pour \alpha = ', num2str(alpha)]);
set(gcf, 'Position', [100, 100, 800, 600]); % Redimensionner la fenêtre de la figure