function y=fbspwavf1(t,m,Fb,Fc)
% fbspwavf morlet
y = (Fb^0.5)*((sinc(Fb*t/m).^m).*exp(2*1i*pi*Fc*t));


%-----------------------------------------------
function y = sinc(x)
%
%               | sin(pi*x)/(pi*x)  if x ~= 0
% y = sinc(x) = |
%               | 1                 if x == 0

y = ones(size(x));
k = find(x);
y(k) = sin(pi*x(k))./(pi*x(k));