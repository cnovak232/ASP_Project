function [xdwt,bands,band_info] = my_dwt(xf,wname,levels,fs,mode,plotting)
% Multi level 1D Discrete Wavelet Transform using matlabs
% 1D 1 level dwt and my own processing 
% computes the frequency at each subband
% computes statistics for each subband 

N = length(xf);
xdwt = zeros(1,N);
bands = struct;
[LoD,HiD,LoR,HiR] = wfilters(wname); 
h0_len = length(LoD);
h1_len = length(HiD);
ext=fix(max(h0_len,h1_len)/2);
xd = xf;

for i = 1:levels
    [xl,xh] = dwt(xd,wname,'mode',mode);
%     xl = xl(1:end-ext+1);
%     xh = xh(1:end-ext+1);
    n = length(xl) + length(xh);
    levelname = ['d' num2str(i)];
    bands.(levelname) = xh; % store level of details
    xdwt(1:n) = [xl,xh];
    xd = xl; 
end
name = 'cA';
bands.(name) = xl; % lowest level band ( approx )

fn = fieldnames(bands);
band_info = struct;
if plotting
    figure;
end
fs_b = fs;
for i=1:length(fn)
    xw = bands.(fn{i});
    if plotting
        subplot(length(fn),1,i);
        plot(xw);
        title(fn{i});
    end
    fs_b = fs_b / 2;
    band_info(i).fs = fs_b;
    band_info(i).means = mean(abs(xh));
    band_info(i).stds = std(abs(xh));
end

end

