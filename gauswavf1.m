function y=gauswavf1(t)
% wavelet support.
% [psi,x] = gauswavf(lb,ub,n,8);
    X2 = t.^2;
    y = (2/pi)^(1/4)*exp(-X2);