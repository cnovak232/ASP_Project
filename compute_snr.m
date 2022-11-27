function snr_db = compute_snr(x,xn)
% x - clean signal
% xn - signal plus an aomunt of noise

resid = x - xn;
snr = mean(x.^2) / mean(resid.^2); % power
snr_db = 10.0 * log10( snr );


end