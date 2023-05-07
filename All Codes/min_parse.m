function [X_r, err] =  min_parse(E, A)

% E = y;
[M, N]= size(A); %M x 5N

if length(E)~=M
    fprintf("Error! Sampled E not correct");
end

%% Normalize A
A_op = A;
for i = 1:N
    A_op(:,i) = A_op(:, i)/norm(A_op(:, i)) ;
end

%% Iteration
s=zeros(1, N);

B = zeros(M, N);

u=zeros(1, N);
r = E;
iteration_size = 500;
err=zeros(1, iteration_size);

for iter = 1:iteration_size
    
    u = abs(A_op'*conj(r));
    
    lambda = find(u == max(u));
    
    test_A = sum(abs(A_op(:, lambda)));
    
    if test_A ~= 0
        s(lambda)=lambda;
        B(:, lambda) = A(:, lambda);
        A_op(:, lambda) = 0;
    end
    
    X_r = (pinv(B' * B)*B')*E;
    r = E - B*X_r;
    
    err_percent = abs(r)./abs(E);
    err(iter) = max((abs(err_percent)));
    
    fprintf("iteration = %d, error = %d, lamdba = %d \n" , iter, err(iter), lambda);
    
    if err(iter)<=1e-5
        break;
    end
end

if length(X_r)~=N
    fprintf("Error! FT not reconstructed properly");
end

end