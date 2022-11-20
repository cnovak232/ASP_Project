function xc = perform_lms(x,r,mu,p)
% LMS Adaptive Noise Cancelation Algorithm
% Inputs
    % x - noisy signal
    % r - reference noise signal
    % mu - step size for convergence
    % p - filter order
% Output
    % xc - noise canceled signal

w = zeros(1,p);
r_buffer = zeros(1,p);
N = length(x);
xc = zeros(N,1);

for n = 1:N
    r_buffer(1) = r(n);
    
    y = w * r_buffer';
    e = x(n) - y;
    xc(n) = e;
    w = w + 2*mu*e.*r_buffer;
    
    r_buffer(end) = 0;
    r_buffer = circshift(r_buffer,1);
end

