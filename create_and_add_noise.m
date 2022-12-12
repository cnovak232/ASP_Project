function [xn,rn] = create_and_add_noise(x,gain,c_order,c_w,type)
% Inputs
% x - primary clean source signal
% gain - level applied to noise
% c_order - the order to filter to simulate channel differences
% type - noise type - gwhite or crowd 

xlen = length(x);
if strcmp(type,'gwhite')
    ref_noise = gain * randn(xlen,1); % noise
else
    [y,fs] = audioread('CrowdNoise/crowd_noise.wav',[1 xlen]);
    y = mean(y,2);
    ref_noise = gain * y;
end

% form some correlation between primary and reference noise
lp = fir1(c_order,c_w);
prim_noise = filter(lp,1,ref_noise); 
xn = x + prim_noise;
rn = ref_noise;

end
