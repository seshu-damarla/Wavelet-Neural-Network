function y=Mexicanhat(t)
% wavelet support.
F = t.^2;
y = (2/(sqrt(3)*pi^0.25)) * exp(-F/2) .* (1-F);
