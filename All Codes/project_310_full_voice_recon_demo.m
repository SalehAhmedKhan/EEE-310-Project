clear all;close all;clc;

%% Voice Signal

[samp_sig_main, Fc] = audioread('myvoice.wav');

samp_sig_main = lowpass(samp_sig_main,3000,Fc);


sound(samp_sig_main, Fc)

sub_sig_main = [];

recon_sig_main = zeros(size(samp_sig_main));

jump = 0;
frame = 1000;
iters = round(length(samp_sig_main)/frame);

%% Frame by Frame Reconstruction

for iter = 1:round(length(samp_sig_main)/frame)
    
    if (frame+jump) >= size(samp_sig_main,1)
        samp_sig = samp_sig_main(1+jump:end, 1);
    else
        samp_sig = samp_sig_main(1+jump:frame+jump, 1);
    end
    
    %% Variables and Constants
    
    N = length(samp_sig) - 1;
    
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
    
    sub_x = mat_main*samp_sig;
    
    sub_sig_main = [sub_sig_main sub_x'];
    
    %% %% Compressive Sensing and Signal Reconstruction
  
    num = 2;
    L=length(samp_sig);
    psy = conj(dftmtx(num*L))/length(samp_sig);
    psy=psy(1:L, :);

    sensing_mat = (mat_main)*psy;
    
    x = fft(samp_sig);
    %figure(3); subplot(3,1,1); stem(fftshift(abs(x))); title('Original signal')
    
    A = sensing_mat;
    y = sub_x;
    %subplot(3,1,2); stem(y); title('Sub sampled values');
    
    % Perform Compressed Sensing recovery
    x0 = A.'*y;
    Xp = l1eq_pd(x0, A, [], y);
    %subplot(3,1,3); stem(fftshift(abs(xp))); title('Recovered signal');
    
    xp = psy * Xp;
    
    %recon_sig = real(xp');
    recon_sig = lowpass(real(xp'),3000,Fc);

    
    %% Plotting
    
    %figure(10)
    %hold on
    %plot(samp_sig,'b.-')
    %plot(circshift(fliplr(real(ifft(xp'))),1),'r.-')
    %legend('nyquist sampled','subsampled reconstructed')
    
    recon_sig_main(1+jump:frame+jump, 1) = recon_sig;
    recon_sig_main(1+jump:frame+jump, 1) = recon_sig;
    
    jump = jump + frame;
    
    %Show Progress
    disp('Iteration: ')
    disp(iter)
    disp(' Out of: ')
    disp(iters)
end

%% Play Audio and Saving
    
%sound(samp_sig_main, Fc)
sound(recon_sig_main, Fc)  

%audiowrite('reconstructed clip 3.wav',recon_sig_main,Fc)

%% SNR Calculation

noise = abs(recon_sig_main(:,1) - samp_sig_main(:,1));

SNR = (sum(samp_sig_main(:,1).*samp_sig_main(:,1))) / (sum(noise.*noise))

SNR_dB = 10*log10(SNR)

%% Plotting

figure(100);
hold on;

super_y = fft(samp_sig_main);
super_x = -Fc/2:Fc/length(samp_sig_main):Fc/2 - Fc/length(samp_sig_main);
stem(super_x, fftshift(abs(super_y(:,1))),'b.');

super_y = fft(recon_sig_main);
super_x = -Fc/2:Fc/length(recon_sig_main):Fc/2 - Fc/length(recon_sig_main);
stem(super_x, fftshift(abs(super_y(:,1))),'r.');

legend('input signal','reconstructed')
title('Signal Comparison (Frequency Domain)')

figure(102)
hold on
plot(samp_sig_main(:,1),'b')
plot(recon_sig_main(:,1),'r')
legend('input signal','reconstructed')
title('Signal Comparison')



    
