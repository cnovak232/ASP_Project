function [xc] = perform_rls(x,r,mu,sigma,lambda,p)
% RLS Adaptive Noise Cancelation Algorithm-
% Sources:
% https://ieeexplore.ieee.org/abstract/document/6016426

% Inputs
    % x - noisy signal
    % r - reference noise signal
    % c - normalization constant (<1)
    % p - filter order
% Output
    % xc - noise canceled signal


w = zeros(1,p);
r_buffer = zeros(1,p);
N = length(x);
xc = zeros(N,1);
P_prev = sigma*eye(p);

for n = 1:N
    r_buffer(1) = r(n);
    
    y = w * r_buffer';
    e = x(n) - y;
    xc(n) = e;
    %calculate weights
    pi = (r_buffer*P_prev) / (lambda + (r_buffer*P_prev*r_buffer'));
    
    w = w + (pi*e);
    P_prev = (lambda^-1)*P_prev - ((lambda^-1)*((pi*r_buffer')*P_prev));

    %J = (lambda^(N-n))*(abs(e)^2);
    
    r_buffer(end) = 0;
    r_buffer = circshift(r_buffer,1);
end