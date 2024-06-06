%% Image loading
x = imread('cameraman.tif');
%x = imresize(x);
[N,M] = size(x);
imshow(x,[0 255])
%% Filter definition
h1 = 1/(25)*ones(5);
H1 = fft2(h1,N,N);

%% creation of Y
SNR = 30;
ynb = conv2(x,h1,"same");
y = adgnoise(ynb,30);
imshow(y,[0 255])

%% Computation of the FFT
Y = fft2(y,N,N);
X_rec = Y./H1;
x_rec = real(ifft2(X_rec));
imshow(x_rec,[0 255]);

%% Use of regularization
alpha = 1;
d1 = [0,-1,0;-1,4,-1;0,-1,0];
D1 = fft2(d1,N,N);
g_mcr = conj(H1)./(abs(H1).^2+alpha*abs(D1).^2);
X_rec2 = g_mcr.*Y;
x_rec2 = real(ifft2(X_rec2));
imshow(x_rec2,[0 255])

%% Choice of alpha

min_alpha=-7;
pas_alpha=0.1;
max_alpha=+3;
i_alpha=0;

x_crop = x(75:175,75:175);
for var_alpha=min_alpha:pas_alpha:max_alpha

alpha=10^var_alpha;
i_alpha=i_alpha+1;
%--------------------------------------------------------------------------
gMCR = conj(H1)./(abs(H1).^2+alpha*abs(D1).^2);
X_hat_regularise = Y .* gMCR;
x_regularise = real(ifft2(X_hat_regularise));
x_regularise = x_regularise(75:175,75:175);
x_rec_l2(:,i_alpha) = x_regularise(:);
%--------------------------------------------------------------------------

diff_err = double(x_crop(:)) - x_rec_l2(:,i_alpha);
err_rec(:,i_alpha)= diff_err;

Werr_rec(i_alpha)= err_rec(:,i_alpha)'*err_rec(:,i_alpha);
end

figure
var_alpha=min_alpha:pas_alpha:max_alpha;
plot(var_alpha,10*log10(Werr_rec))
xlabel('log10(\alpha)'); ylabel('energy of reconstruction error (dB)'); title('Reconstruction error with respect to \alpha');

%% y
figure(2)
imshowpair(x,y,"montage") 


%% Optimization
clear options

alpha = 0.4;
T = 4;
x0 = y;
f = @(x)F_alpha(x,y,h1,alpha,T);
options = optimoptions('fminunc', 'HessianApproximation','lbfgs','SpecifyObjectiveGradient',true, 'OptimalityTolerance',0, 'StepTolerance',0, 'FunctionTolerance', 0, 'MaxIterations',10);
options.CheckGradients = false;
x_rec_1 = fminunc(f,x0,options);
diff_err_1 = double(x(:)) - x_rec_1(:);
err_rec_1 = diff_err_1'*diff_err_1/(256*256);
figure(2);
imshow(x_rec_1, [0, 255])
%%
figure
tiledlayout(1,5)

%nexttile
%imshow(y, [0, 255])

alpha = 1;
T = 10;
x0 = y;
f = @(x)F_alpha(x,y,h1,alpha,T);
options = optimoptions('fminunc', 'HessianApproximation','lbfgs','SpecifyObjectiveGradient',true, 'OptimalityTolerance',0, 'StepTolerance',0, 'FunctionTolerance', 0, 'MaxIterations',10);
options.CheckGradients = false;
x_rec_1 = fminunc(f,x0,options);
diff_err_1 = double(x(:)) - x_rec_1(:);
err_rec_1 = diff_err_1'*diff_err_1/(256*256);
nexttile
imshow(x_rec_1, [0, 255])
title(err_rec_1)

alpha = 1;
T = 10;
x0 = y;
f = @(x)F_alpha(x,y,h1,alpha,T);
options = optimoptions('fminunc', 'HessianApproximation','lbfgs','SpecifyObjectiveGradient',true, 'OptimalityTolerance', 0, 'StepTolerance',0, 'FunctionTolerance', 0, 'MaxIterations',20);
options.CheckGradients = false;
x_rec_2 = fminunc(f,x0,options);
diff_err_2 = double(x(:)) - x_rec_2(:);
err_rec_2 = diff_err_2'*diff_err_2/(256*256);
nexttile
imshow(x_rec_2, [0, 255])
title(err_rec_2)

alpha = 1;
T = 10;
x0 = y;
f = @(x)F_alpha(x,y,h1,alpha,T);
options = optimoptions('fminunc', 'HessianApproximation','lbfgs','SpecifyObjectiveGradient',true,'OptimalityTolerance',0, 'StepTolerance',0, 'FunctionTolerance', 0, 'MaxIterations',40);
options.CheckGradients = false;
x_rec_3 = fminunc(f,x0,options);
diff_err_3 = double(x(:)) - x_rec_3(:);
err_rec_3 = diff_err_3'*diff_err_3/(256*256);
nexttile
imshow(x_rec_3, [0, 255])
title(err_rec_3)


alpha = 1;
T = 10;
x0 = y;
f = @(x)F_alpha(x,y,h1,alpha,T);
options = optimoptions('fminunc', 'HessianApproximation','lbfgs','SpecifyObjectiveGradient',true, 'OptimalityTolerance', 0, 'StepTolerance',0, 'FunctionTolerance', 0, 'MaxIterations',80);
options.CheckGradients = false;
x_rec_4 = fminunc(f,x0,options);
diff_err_4 = double(x(:)) - x_rec_4(:);
err_rec_4 = diff_err_4'*diff_err_4/(256*256);
nexttile
imshow(x_rec_4, [0, 255])
title(err_rec_4)

alpha = 1;
T = 10;
x0 = y;
f = @(x)F_alpha(x,y,h1,alpha,T);
options = optimoptions('fminunc', 'HessianApproximation','lbfgs','SpecifyObjectiveGradient',true, 'OptimalityTolerance',0, 'StepTolerance',0, 'FunctionTolerance', 0, 'MaxIterations',2000);
options.CheckGradients = false;
x_rec_5 = fminunc(f,x0,options);
diff_err_5 = double(x(:)) - x_rec_5(:);
err_rec_5 = diff_err_5'*diff_err_5/(256*256);
nexttile
imshow(x_rec_5, [0, 255])
title(err_rec_5)


%figure(2);
%montage([y,x_rec_1, x_rec_2, x_rec_3, x_rec_4], "Displayrange", [0, 255]) 

function [err, grad] = F_alpha(x, y, h, alpha, T)
    [M, N] = size(x);
    conv_result = conv2(x, h, "same");
    term1 = sum((y(:) - conv_result(:)).^2);
    
    term2 = 0;
    grad_term = zeros(M, N);
    
    for i = 2:M-1
        for j = 2:N-1
            % Compute differences
            diff1 = x(i,j) - x(i-1,j);
            diff2 = x(i,j) - x(i+1,j);
            diff3 = x(i,j) - x(i,j-1);
            diff4 = x(i,j) - x(i,j+1);
            
            % Compute phi value
            phi = sqrt(diff1^2 + T^2) + sqrt(diff2^2 + T^2) + sqrt(diff3^2 + T^2) + sqrt(diff4^2 + T^2) - 4*T;
            term2 = term2 + phi;
            
            % Gradient contributions from regularization term
            grad_term(i,j) = grad_term(i,j) + diff1/sqrt(diff1^2 + T^2) + diff2/sqrt(diff2^2 + T^2) + diff3/sqrt(diff3^2 + T^2) + diff4/sqrt(diff4^2 + T^2);
            
            % Adjust the neighboring pixels' gradients accordingly
            grad_term(i-1,j) = grad_term(i-1,j) - diff1/sqrt(diff1^2 + T^2);
            grad_term(i+1,j) = grad_term(i+1,j) - diff2/sqrt(diff2^2 + T^2);
            grad_term(i,j-1) = grad_term(i,j-1) - diff3/sqrt(diff3^2 + T^2);
            grad_term(i,j+1) = grad_term(i,j+1) - diff4/sqrt(diff4^2 + T^2);
        end
    end
    
    term2 = alpha * term2;
    err = term1 + term2;
    
    % Gradient of data fidelity term
    grad_data_fidelity = -2 * conv2(y - conv_result, h', "same");
    
    % Combine gradients
    grad = grad_data_fidelity + alpha * grad_term;
end