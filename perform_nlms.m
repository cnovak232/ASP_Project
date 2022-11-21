function xc = perform_nlms(x,r,alpha,p)
% LMS Adaptive Noise Cancelation Algorithm
% Inputs
    % x - noisy signal
    % r - reference noise signal
    % alpha - adaption constant 0<alpha<2
    % p - filter order
% Output
    % xc - noise canceled signal

w = zeros(1,p);
r_buffer = zeros(1,p);
N = length(x);
xc = zeros(N,1);
c = .001; % safety factor

for n = 1:N
    r_buffer(1) = r(n);
    
    y = w * r_buffer';
    e = x(n) - y;
    xc(n) = e;
    mu =  2 * alpha / (c + norm(r_buffer));
    w = w + mu*e.*r_buffer;
    
    r_buffer(end) = 0;
    r_buffer = circshift(r_buffer,1);
end