function psnr_db = compute_psnr(x,xn)

resid = x - xn;
psnr = max(x).^2 / mean(resid.^2); % power
psnr_db = 10.0 * log10( psnr );

end

