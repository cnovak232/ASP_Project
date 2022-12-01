%% ASP Term Project Driver
% Adaptive Noise Cancelation in Music signals
% Chris Novak, Natalie Meyer

%% Read in audio data
audiodir = './ASP_Project_Audio/';
listname = dir(audiodir);
listname = listname(3:end);
fs = 44100;
t_per_song = 5; % 10 second clips of each song
num_samples = t_per_song * fs;
music_files = {};
for i = 1:length(listname)
    [y,fs] = audioread([audiodir, listname(i).name],[1 num_samples]);
    music_files{i} = y;
end


%% Add variable noise noise
x = music_files{1};

% make mono for now
x = mean(x,2);

xlen = length(x);
noiselen = 1 * fs; % change noise every 1 sec
num_frames = xlen / noiselen;
xn = zeros(xlen,1); % signal + noise
ref_noise = zeros(xlen,1); % noise
gains = [.08 .08 .08 .08 .08]; % can change this to be varying
for f = 0:num_frames-1
    ind = 1+f*noiselen:noiselen+f*noiselen;
    ref_noise(ind) = gains(f+1)*randn(noiselen,1);
end

% form some correlation between primary and reference noise
lp = fir1(11,.4);
prim_noise = filter(lp,1,ref_noise); 
xn = x + prim_noise;

figure;
subplot(311);
plot(x);
title('Clean Signal');
subplot(312);
plot(ref_noise); 
title('Reference Noise Signal')
subplot(313);
plot(xn);
title('Signal with noise')

%% Run adaptive noise cancelation algorithms

xc_lms = perform_lms(xn,ref_noise,.001,10);
xc_nlms = perform_nlms(xn,ref_noise,.001,10);
xc_rls = perform_rls(xn,ref_noise,1,1,10);
xc_afa = perform_afa(xn,ref_noise,.6,10);

% figure;
% subplot(411);
% plot(x);
% subplot(412);
% plot(xn);
% subplot(413);
% plot(xc_lms);
% subplot(414);
% plot(xc_nlms);

% plot convergence of algorithms
converge_lms = abs(x - xc_lms);
converge_nlms = abs(x - xc_nlms);
converge_rls = abs( x - xc_rls);
converge_afa = abs(x - xc_afa);
figure;
subplot(411);
plot(converge_lms);
subplot(412);
plot(converge_nlms);
subplot(413);
plot(converge_rls);
subplot(414);
plot(converge_afa);

% Compare SNR 
snr_before = compute_snr(x,xn)

snr_lms = compute_snr(x,xc_lms)

snr_nlms = compute_snr(x,xc_nlms)

snr_rls = compute_snr(x,xc_rls)

snr_afa = compute_snr(x,xc_afa)

%% Multirate

xc_wavelet = perform_wavelet_anc(xn,ref_noise,.01,10,1,'db6');

snr_mr = compute_snr(x,xc_wavelet)
