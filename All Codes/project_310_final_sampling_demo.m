clear all;close all;clc;

%% Variables and Constants
Fc = 1000;

F1 = 91;
F2 = 103;
F3 = 193;

N1 = Fc/F1;
N2 = Fc/F2;
N3 = Fc/F3;
N =  1000 - 1;
M1 = floor(F1*(N/Fc));
M2 = floor(F2*(N/Fc))-1;
M3 = floor(F3*(N/Fc))-1;
M = M1+M2+M3;

delay12 = 1/F2;
delay13 = 1/F3;

Tc = 1/Fc;

%% Signals
sig = @(t) sin(2*pi*70*t) + sin(2*pi*275*t) + sin(2*pi*400*t);
%sig = @(t) sin(2*pi*100*t);
%sig = @(t) sinc(30*t - 15);
%sig = @(t) rectangularPulse(0.35,0.65,t);

%% Sampling
figure(1)
hold on

%ANALOG SIGNAL
start_t = 0;
end_t = N*Tc;
t = start_t:0.0001:end_t;
analog_sig = sig(t);
plot(t,analog_sig);
%plot(t,analog_sig,'--');

%SAMPLING (NYQUIST RATE = 800, SAMPLING RATE = 1000)
Nc = 0:N;
Nc_time = Nc/Fc;
samp_sig = sig(Nc_time);
%recon_sig = ws_interp(samp_sig, Fc, t);
%plot(t,recon_sig,'k-.');

%SUB NYQUIST SAMPLING
Ns1 = 0:M1;
Ns1_time = Ns1/F1;
sub_sig1 = sig(Ns1_time);

Ns2 = 0:M2;
Ns2_time = Ns2/F2  + delay12;
sub_sig2 = sig(Ns2_time);

Ns3 = 0:M3;
Ns3_time = Ns3/F3 + delay13;
sub_sig3 = sig(Ns3_time);

% stem(Nc_time,samp_sig);
% stem(Ns1_time,sub_sig1,'r');
% stem(Ns2_time,sub_sig2,'b');
% stem(Ns3_time,sub_sig3,'k');

%stem(Nc_time,samp_sig);
stem(Ns1_time,sub_sig1,'r*-');
stem(Ns2_time,sub_sig2,'b*-');
stem(Ns3_time,sub_sig3,'k*-');
%legend('analog','nyquist sample')
title('Sampled Signals')


%% Measurement Matrices
mat1 = zeros(M1, N);
mat2 = zeros(M2, N);
mat3 = zeros(M3, N);

for m = 0:M1
    for n = 0:N
        mat1(m+1,n+1) = sinc((m*N1*Tc)/Tc - n);
    end
end

for m = 0:M2
    for n = 0:N
        mat2(m+1,n+1) = sinc((m*N2*Tc + delay12)/Tc - n);
    end
end

for m = 0:M3
    for n = 0:N
        mat3(m+1,n+1) = sinc((m*N3*Tc + delay13)/Tc - n);
    end
end

mat_main = [mat1' mat2' mat3']';

sub_x2 = mat_main*samp_sig';
sub_x = [sub_sig1 sub_sig2 sub_sig3]';


% figure(2)
% hold on
% stem(sub_x,'b');
% stem(sub_x2,'r*');
% legend('theory','actual');

%% Compressive Sensing and Signal Reconstruction

num = 1;
L=length(samp_sig);

psy = conj(dftmtx(num*L))/length(samp_sig);
psy=psy(1:L, :);

sensing_mat = (mat_main)*psy;

x = (fft(samp_sig, num*L));
figure(3); 
subplot(2,1,1); 
stem(fftshift(abs(x)),'.'); 
title('Original Fourier Spectra')

%% Time Sorting

time_axis=[Ns1_time, Ns3_time, Ns2_time]';
sort_mat = zeros (M+3, 2 + L);
sort_mat(:, 1) = time_axis;
sort_mat(:, 2) = sub_x;
sort_mat(:, 3:end) = sensing_mat;

sort_mat = sortrows(sort_mat,1);

sub_x = sort_mat(:, 2);
sensing_mat = sort_mat(:, 3:end);


A = sensing_mat;
y = sub_x;
%y = sub_x2;
%subplot(3,1,2);stem(y,'.');title('Sub sampled values'); 

% Perform Compressed Sensing Recovery
x0 = A.'*y;
Xp = l1eq_pd(x0, A, [], y);

subplot(2,1,2); stem(fftshift(abs(Xp)),'.'); title('Recovered FT');

%% Results

figure(15)
hold on
plot(t,analog_sig,'b');
recon_sig = ws_interp(samp_sig, Fc, t);
%plot(t,recon_sig,'b');

xp = psy * Xp;
recon_sig_2 = ws_interp(real(xp'), Fc, t);
plot(t,recon_sig_2,'r');
legend('analog signal','subsampled reconstructed')
title('Signal Comparison');

figure(16)
hold on
stem(samp_sig,'b.')
stem(real(xp'),'r.')
legend('input signal','subsampled reconstructed')
title('Signal Comparison (Discrete)');


%% SNR

recon_sig_main = real(xp');
%recon_sig_main = recon_sig_2;

samp_sig_main = samp_sig;
%samp_sig_main = analog_sig;

noise = abs(recon_sig_main - samp_sig_main);

SNR = (sum(samp_sig_main.*samp_sig_main)) / (sum(noise.*noise))

SNR_dB = 10*log10(SNR)











