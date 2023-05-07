clear all;close all;clc;

%% Voice Signal
[samp_sig_main, Fc] = audioread('myvoice.wav');

samp_sig_main = lowpass(samp_sig_main,3000,Fc); %Bandlimiting the Input Signal

jump = 43000;

samp_sig = samp_sig_main(1+jump:1000+jump);

%% Variables and Constants

N = size(samp_sig, 2) - 1;

Tc = 1/Fc;

F1 = 3300;
F2 = 4650;
F3 = 7430;

N1 = Fc/F1;
N2 = Fc/F2;
N3 = Fc/F3;
M1 = floor(F1*(N/Fc));
M2 = floor(F2*(N/Fc));
M3 = floor(F3*(N/Fc));
M = M1+M2+M3;

delay12 = 1/F3;
delay13 = 1/F2;

%% Measurement Matrices

mat1 = zeros(M1, N);
mat2 = zeros(M2, N);
mat3 = zeros(M3, N);

for m = 0:M1
    for n = 0:N
        mat1(m+1,n+1) = sinc((m*N1) - n);
    end
end

for m = 0:M2
    for n = 0:N
        mat2(m+1,n+1) = sinc((m*N2*Tc - delay12)/Tc - n);
    end
end

for m = 0:M3
    for n = 0:N
        mat3(m+1,n+1) = sinc((m*N3*Tc - delay13)/Tc - n);
    end
end

mat_main = [mat1' mat2' mat3']';

sub_x = mat_main*samp_sig';

%% Compressive Sensing and Signal Reconstruction

num = 5;
L=length(samp_sig);

psy = conj(dftmtx(num*L))/length(samp_sig);
psy=psy(1:L, :);

sensing_mat = (mat_main)*psy;

x = fft(samp_sig, L*num);
figure(3); subplot(2,1,1); stem(fftshift(abs(x)),'.'); title('Original Fourier Spectra')

A = sensing_mat;
y = sub_x;
%subplot(3,1,2); stem(y,'.'); title('Sub sampled values'); 

% Perform Compressed Sensing recovery
x0 = A.'*y;
Xp = l1eq_pd(x0, A, [], y);
subplot(2,1,2); stem(fftshift(abs(Xp')),'.'); title('Recovered FT');

xp = psy * Xp;

%recon_sig_main = real(xp');
recon_sig_main = lowpass(real(xp'),3000,Fc); %Bandlimiting the Reconstructed Signal

%% Plotting

figure(10)
hold on
plot(samp_sig,'b.-')
plot(recon_sig_main,'r.-')
legend('nyquist sampled','subsampled reconstructed')
title('Signal Comparison')

%% SNR Calculation

recon_sig_main = real(xp');
samp_sig_main = samp_sig;

noise = abs(recon_sig_main - samp_sig_main);

SNR = (sum(samp_sig_main.*samp_sig_main)) / (sum(noise.*noise));

SNR_dB = 10*log10(SNR)



















