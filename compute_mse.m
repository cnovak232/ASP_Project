function mse = compute_mse(x,xc)
% x - desired signal
% cleaned/reconstructed signal

mse = mean( (x-xc).^2 );

end

