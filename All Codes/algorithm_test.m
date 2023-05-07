close all; clear all; clc;

% Initialize constants and variables
rng(0);                 % set RNG seed
N = 256;                % length of signal
P = 40;                  % number of non-zero peaks
K = 100;                 % number of measurements to take (N < L)
x = zeros(N,1);         % original signal (P-sparse)

% Generate signal with P randomly spread values
peaks = randperm(N);
peaks = peaks(1:P);
x(peaks) = randn(1, P);
amp = 1.2*max(abs(x));
figure; subplot(3,1,1); plot(x); title('Original signal'); xlim([1 N]); ylim([-amp amp]);

% Obtain K measurements
A = randn(K, N);
y = A*x;
subplot(3,1,2); plot(y); title('K measured values'); xlim([1 K]);

%% Perform Compressed Sensing recovery
% x0 = A'*y;
% xp = l1eq_pd(x0, A, [], y);
% subplot(3,1,3); plot(real(xp)); title('Recovered signal'); xlim([1 N]); ylim([-amp amp]);

%% Sparse Algorithm

[X_r, err] = min_parse(y, A);

% E = y;
% [M, len_5N]= size(A); %M x 5N
% 
% if length(E)~=M
%     fprintf("Error! Sampled E not correct");
% end
% 
% %% Normalize A
% 
% A_op = A;
% for i = 1:len_5N
%     A_op(:,i) = A_op(:, i)/norm(A_op(:, i)) ;
% end
% 
% %% Iteration
% s=zeros(1, len_5N);
% 
% B = zeros(M, len_5N);
% % B = randn(M, len_5N);
% % B = (B/max(max(abs(B))))*1e-5;
% 
% u=zeros(1, len_5N);
% r = E;
% iteration_size = 500;
% err=zeros(1, iteration_size);
% 
% for iter = 1:iteration_size
%     
%     u = abs(A_op'*r);
%     
%     lambda = find(u == max(u));
%     
%     test_A = sum(abs(A_op(:, lambda)));
%     
%     if test_A ~= 0
%         s(lambda)=lambda;
%         B(:, lambda) = A(:, lambda);
%         A_op(:, lambda) = 0;
%     end
%     
%     X_r = (pinv(B' * B)*B')*E;
%     r = E - B*X_r;
%     
%     err_percent = abs(r)./abs(E);
%     err(iter) = max((abs(err_percent)));
%     
%     fprintf("iteration = %d, error = %d, lamdba = %d \n" , iter, err(iter), lambda);
%     
%     if err(iter)<=1e-5
%         break;
%     end
% end
% 
% % err = real(err);
% 
% if length(X_r)~=len_5N
%     fprintf("Error! FT not reconstructed properly");
% end

%Plot
xp = X_r;
subplot(3,1,3); plot(real(xp)); title('Recovered signal'); xlim([1 N]); ylim([-amp amp]);

figure(2);
plot(err);

err_recons= xp-x;
SNR_dB = 10*log10((x'*x)/(err_recons'*err_recons))

