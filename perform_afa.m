function [xc] = perform_afa(x,r,gamma,p)
% Adaptive Average Filter Noise Cancelation Algorithm
%Sources:
%https://ieeexplore-ieee-org.proxy1.library.jhu.edu/document/9418193
%https://ieeexplore.ieee.org/document/1300608
% Inputs
    % x - noisy signal
    % r - reference noise signal
    % p - filter order
    % gamma - gain parameter
% Output
    % xc - noise canceled signal

w = zeros(p,1);
r_buffer = zeros(p,1);
N = length(x);
xc = zeros(N,1);
wprev = zeros(p,1);
gprev = zeros(p,1);

for n = 1:N
    r_buffer(1) = r(n);
    
    y = w' * r_buffer;
    e = x(n) - y;
    xc(n) = e;
   % w = w + 2*mu*e.*r_buffer;

    what = (1/n) .* (((n-1).*wprev) + w);
    ghat = (1/(n^gamma)) * ((((n-1)^gamma) * gprev) + e.*r_buffer);

    w = what + ghat;

    wprev = what;
    gprev = ghat;
    
    r_buffer(end) = 0;
    r_buffer = circshift(r_buffer,1);
end
end

