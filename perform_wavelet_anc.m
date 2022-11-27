function xc = perform_wavelet_anc(x,r,mu,p,lvl,wname)
% This version of multirate adaptive processing applies wavelet based
% thresholding to the detail bands (number depends on level of wavelet iters) 
% while performing adaptive lms on the lowest approximation band. It then
% reconstructions the signal for a cleaner signal 

N = length(x);
xdwt = zeros(1,N);
rdwt = zeros(1,N);
bands_x = struct;
bands_r = struct;

xd = x;
rd = r;

% perform multilevel wavelet decomposition
for i = 1:lvl
    [xl,xh] = dwt(xd,wname,'mode','per');
    [rl,rh] = dwt(rd,wname,'mode','per');
    nx = length(xl) + length(xh);
    nr = length(rl) + length(rh);
    levelname = ['d' num2str(i)];
    bands_x.(levelname) = xh; % store level of details
    bands_r.(levelname) = rh;
    xdwt(1:nx) = [xl,xh];
    rdwt(1:nr) = [rl,rh];
    xd = xl;
%     subplot(lvl+1,1,i);
%     plot(xh);
%     title(levelname);
end

% run adaptive lms on the lowest band (coarse appriximation)
xc_l = perform_lms(xl,rl,mu,p);

fn = fieldnames(bands_x);

% threshold the detail bands (ideally  mostly noise).
for i = 1:length(fn)
    band = bands_x.(fn{i});

    thr = thselect(band,'rigrsure');
    
    band_thr = wthresh(band,'s',thr);
    
    bands_x.(fn{i}) = band_thr;
end

% perfrom multi level idwt to reconstruction signal
xD = bands_x.(fn{end});
xA = xc_l;
for i = 1:lvl
    xc_mr = idwt(xA,xD,wname,'mode','per');
    if (i~=lvl)
        xA = xc_mr;
        xD = bands_x.(fn{end-i});
    end
end

xc = xc_mr;