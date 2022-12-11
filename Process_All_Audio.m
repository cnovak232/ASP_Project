%% Script for looping through all audio files and processing

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
%% Loop through audio clips and compute average PSNR improvement
psnr_lms = zeros(1,length(p));
psnr_nlms = zeros(1,length(p));
psnr_rls = zeros(1,length(p));
psnr_afa = zeros(1,length(p));
for i = 1:length(music_files)
x = music_files{i};

% make mono for now
x = mean(x,2);

xlen = length(x);
noiselen = 1 * fs; % change noise every 1 sec
num_frames = xlen / noiselen;
ref_noise = zeros(xlen,1); % noise
gains = [.1 .1 .1 .1 .1 .1 .1 .1 .1 .1]; % can change this to be varying
for f = 0:num_frames-1
    ind = 1+f*noiselen:noiselen+f*noiselen;
    ref_noise(ind) = gains(f+1)*randn(noiselen,1);
end

% form some correlation between primary and reference noise
lp = fir1(10,.4);
prim_noise = filter(lp,1,ref_noise); 
xn = x + prim_noise;
p = 10;
xc_lms = perform_lms(xn,ref_noise,best_mu_lms(7),p);
xc_nlms = perform_nlms(xn,ref_noise,best_mu_nlms(7),p);
xc_rls = perform_rls(xn,ref_noise,best_lam_rls(7),1,p);
xc_afa = perform_afa(xn,ref_noise,best_gam_afa(7),p); 

% plot convergence of algorithms
converge_lms = abs(x - xc_lms);
converge_nlms = abs(x - xc_nlms);
converge_rls = abs( x - xc_rls);
converge_afa = abs(x - xc_afa);
% figure;
% subplot(411);
% plot(converge_lms);
% subplot(412);
% plot(converge_nlms);
% subplot(413);
% plot(converge_rls);
% subplot(414);
% plot(converge_afa);

% Compare SNR 
psnr_before = compute_psnr(x,xn);

psnr_lms(i) = compute_psnr(x,xc_lms) - psnr_before;

psnr_nlms(i) = compute_psnr(x,xc_nlms)- psnr_before;

psnr_rls(i) = compute_psnr(x,xc_rls)- psnr_before;

psnr_afa(i) = compute_psnr(x,xc_afa)- psnr_before;

end

psnr_lms_avg = mean(psnr_lms)
psnr_nlms_avg = mean(psnr_nlms)
psnr_rls_avg = mean(psnr_rls)
psnr_afa_avg = mean(psnr_afa)


% audiowrite('clean_x.wav',x,fs);
% audiowrite('unclean_x.wav',xn,fs);
% 
% audiowrite('clean_lms.wav',xc_lms,fs);
% audiowrite('clean_nlms.wav',xc_nlms,fs);
% audiowrite('clean_rls.wav',xc_rls,fs);
% audiowrite('clean_afa.wav',xc_afa,fs); 
