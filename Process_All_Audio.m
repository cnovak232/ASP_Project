%% Script for looping through all audio files and processing

%% Read in audio data
audiodir = './ASP_Project_Audio/';
listname = dir(audiodir);
listname = listname(3:end);
fs = 44100;
t_per_song = 10; % 10 second clips of each song
num_samples = t_per_song * fs;
music_files = {};
for i = 1:length(listname)
    [y,fs] = audioread([audiodir, listname(i).name],[1 num_samples]);
    music_files{i} = y;
end

%% 
audiodir = './SpeechSignals/';
listname = dir(audiodir);
listname = listname(3:end);
 fs = 16000;
 t_per_song = 10; % 10 second clips of each song
 num_samples = t_per_song * fs;
speech_files = {};
for i = 1:length(listname)
    [y,fs] = audioread([audiodir, listname(i).name]);
    speech_files{i} = y;
end

%% Loop through audio clips and compute average PSNR improvement
audio_files = music_files; % music_files
len = length(audio_files);
snr_lms = zeros(1,len);
snr_nlms = zeros(1,len);
snr_rls = zeros(1,len);
snr_afa = zeros(1,len);

for i = 1:len

x = audio_files{i};

% make mono for now
x = mean(x,2);

[xn,ref_noise] = create_and_add_noise(x,.5,10,.5,'crowd');

p = 10;
xc_lms = perform_lms(xn,ref_noise,best_params.mu_lms,p);
xc_nlms = perform_nlms(xn,ref_noise,best_params.mu_nlms,p);
xc_rls = perform_rls(xn,ref_noise,best_params.lam_rls,1,p);
xc_afa = perform_afa(xn,ref_noise,best_params.gam_afa,p); 

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
snr_before = compute_snr(x,xn);

snr_lms(i) = compute_snr(x,xc_lms) - snr_before;

snr_nlms(i) = compute_snr(x,xc_nlms)- snr_before;

snr_rls(i) = compute_snr(x,xc_rls)- snr_before;

snr_afa(i) = compute_snr(x,xc_afa)- snr_before;

end

psnr_lms_avg = mean(snr_lms)
psnr_nlms_avg = mean(snr_nlms)
psnr_rls_avg = mean(snr_rls)
psnr_afa_avg = mean(snr_afa)


% audiowrite('clean_x.wav',x,fs);
% audiowrite('unclean_x.wav',xn,fs);
% 
% audiowrite('clean_lms.wav',xc_lms,fs);
% audiowrite('clean_nlms.wav',xc_nlms,fs);
% audiowrite('clean_rls.wav',xc_rls,fs);
% audiowrite('clean_afa.wav',xc_afa,fs); 
