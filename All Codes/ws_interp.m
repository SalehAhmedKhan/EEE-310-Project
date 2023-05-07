function res = ws_interp(samp,Fs,t)

N = length(samp);
T = 1/Fs;

res = 0;
for i=1:N
    res = res + samp(i) * sinc((t-(i-1)*T)/T);
end
