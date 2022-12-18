function xc = perform_wavelet_anc(x,r,anc,anc_param,p,lvl,wname,fs)
% This version of multirate adaptive processing applies wavelet based
% thresholding to the detail bands (number depends on level of wavelet iters) 
% while performing adaptive lms on the lowest approximation band. It then
% reconstructions the signal for a cleaner signal 
mode = 'per';
[xdwt,xbands,xband_info] = my_dwt(x,wname,lvl,fs,mode,0);
[rdwt,rbands, rband_info] = my_dwt(r,wname,lvl,fs,mode,0);

% run adaptive lms on the lowest band (coarse appriximation)
xc_cA = anc(xbands.cA,rbands.cA,anc_param,p);

fn = fieldnames(xbands);

% soft threshold the detail bands (ideally  mostly noise).
for i = 1:length(fn)-1
    band = xbands.(fn{i});

%     xwt_sort = sort(abs(xwt(:)),'descend');
%     ind = ceil( 0.10 * length(xwt_sort));
%     thresh = xwt_sort(ind);
%     
%     xwt(abs(xwt) < thresh ) = 0;
    thr = thselect(band,'rigrsure');
    
    band_thr = wthresh(band,'s',thr);
    
    xbands.(fn{i}) = band_thr;
end

% perfrom multi level idwt to reconstruction signal
xD = xbands.(fn{end-1});
xA = xc_cA;
for i = 1:lvl
    xc_mr = idwt(xA,xD,wname,'mode',mode);
    if (i~=lvl)
        xD = xbands.(fn{(end-1)-i});
        xA = xc_mr(1:length(xD));
        %xD = [xD; zeros(length(xA)-length(xD),1)];
    end
end

xc = xc_mr;