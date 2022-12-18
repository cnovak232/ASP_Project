function ser = compute_ser(x,xc)
% x - desired signal
% cleaned/reconstructed signal

mse = mean(xc.^2) / mean( (x-xc).^2 );
ser = 10.0* log10(mse);

end

