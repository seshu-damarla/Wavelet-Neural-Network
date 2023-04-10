function y=d_fbspwavf1(t,m,Fb,Fc)
% wavelet support.
% [psi,x] = gauswavf(lb,ub,n,8);
  y= Fb^(1/2)*Fc*pi*exp(pi*Fc*t*2i)*((m*sin((pi*Fb*t)/m))/(Fb*t*pi))^m*2i + Fb^(1/2)*m*exp(pi*Fc*t*2i)*(cos((pi*Fb*t)/m)/t - (m*sin((pi*Fb*t)/m))/(Fb*t^2*pi))*((m*sin((pi*Fb*t)/m))/(Fb*t*pi))^(m - 1);