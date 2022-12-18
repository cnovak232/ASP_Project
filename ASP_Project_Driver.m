%% ASP Term Project Driver
% Adaptive Noise Cancelation in Music signals
% Chris Novak, Natalie Meyer

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
%% Read in Speech Signals 

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

%% Add variable noise noise
x = music_files{1};

% make mono for now
x = mean(x,2);

[xn,ref_noise] = create_and_add_noise(x,.1,10,.5,'gwhite');

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

%% Run adaptive noise cancelation algorithms - optimal params
p = 10; % filter order
mu = .001; % convergence factor for lms/nlms between 0 and 1
lambda = 1; % "forgetting" factor for rls - usually between .98 and 1
gamma = .5; % gain parameter for afa convergence between .5 and 1

xc_lms = perform_lms(xn,ref_noise,.0231,p);
xc_nlms = perform_nlms(xn,ref_noise,.002,p);
xc_rls = perform_rls(xn,ref_noise,1,p);
xc_afa = perform_afa(xn,ref_noise,.5,p);

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
ylabel(' Error');
title('LMS Convergence')
xlabel('Samples n (iterations)')
subplot(412);
plot(converge_nlms);
ylabel('Error');
title('NLMS Convergence')
xlabel('Samples n (iterations)')
subplot(413);
plot(converge_rls);
ylabel('Error');
title('RLS Convergence')
xlabel('Samples n (iterations)')
ylim([0 .2]);
subplot(414);
plot(converge_afa);
ylim([0 .2]);
ylabel('Error');
title('AFA Convergence')
xlabel('Samples n (iterations)')
sgtitle('Convergence of Algorithms with Optimal Parameters')

% Compare MSE
mse_before = compute_ser(x,xn)

mse_lms = compute_ser(x,xc_lms) - mse_before

mse_nlms = compute_ser(x,xc_nlms) - mse_before

mse_rls = compute_ser(x,xc_rls) - mse_before

mse_afa = compute_ser(x,xc_afa) - mse_before


Compare SNR 
snr_before = compute_snr(x,xn)

snr_lms = compute_snr(x,xc_lms)

snr_nlms = compute_snr(x,xc_nlms)

snr_rls = compute_snr(x,xc_rls)

snr_afa = compute_snr(x,xc_afa)

% Compare PSNR 
psnr_before = compute_psnr(x,xn)

psnr_lms_imp = compute_psnr(x,xc_lms) - psnr_before

psnr_nlms_imp = compute_psnr(x,xc_nlms) - psnr_before

psnr_rls_imp = compute_psnr(x,xc_rls) - psnr_before

psnr_afa_imp = compute_psnr(x,xc_afa) - psnr_before

%% Multirate
p = 10;
levels = 1;
wname = 'db8';
xc_dwt_lms = perform_wavelet_anc(xn,ref_noise,@perform_lms,.0231,p,levels,'db6',fs);
xc_dwt_nlms = perform_wavelet_anc(xn,ref_noise,@perform_nlms,.002,p,levels,'db6',fs);
xc_dwt_rls = perform_wavelet_anc(xn,ref_noise,@perform_rls,1,10,levels,'db6',fs);
xc_dwt_afa = perform_wavelet_anc(xn,ref_noise,@perform_afa,.5,10,levels,'db6',fs);

% Compare MSE
mse_before = compute_ser(x,xn)

mse_lms = compute_ser(x,xc_dwt_lms) - mse_before

mse_nlms = compute_ser(x,xc_dwt_nlms) - mse_before

mse_rls = compute_ser(x,xc_dwt_rls) - mse_before

mse_afa = compute_ser(x,xc_dwt_afa) - mse_before


snr_mr = compute_snr(x,xc_dwt_lms)

% xd = decimate(xn,2,'fir');
% rd = decimate(ref_noise,2,'fir');
% 
% xc_d_lms = perform_lms(xd,ref_noise,.0231,10);
% 
% xc_i_lms = interp(xc_d_lms,2);
% 
% mse_d = 10*log10(compute_mse(x,xc_i_lms))
% snr_d = compute_snr(x,xc_i_lms)

%% Compute optimal convergence params for given system

best_params = compute_optimal_params(x,xn,ref_noise,10);


%% Compare Filter order vs SNR

p = 4:32;
mse_lms = zeros(1,length(p));
mse_nlms = zeros(1,length(p));
mse_rls = zeros(1,length(p));
mse_afa = zeros(1,length(p));

snr_lms = zeros(1,length(p));
snr_nlms = zeros(1,length(p));
snr_rls = zeros(1,length(p));
snr_afa = zeros(1,length(p));

psnr_lms = zeros(1,length(p));
psnr_nlms = zeros(1,length(p));
psnr_rls = zeros(1,length(p));
psnr_afa = zeros(1,length(p));

for i = 1:length(p)
    xc_lms = perform_lms(xn,ref_noise,best_mu_lms,p(i));
    xc_nlms = perform_nlms(xn,ref_noise,best_mu_nlms,p(i));
    xc_rls = perform_rls(xn,ref_noise,best_lam_rls,1,p(i));
    xc_afa = perform_afa(xn,ref_noise,best_gam_afa,p(i));

    mse_lms(i) = 10*log10( compute_mse(x,xc_lms) );
    mse_nlms(i) = 10*log10( compute_mse(x,xc_nlms));
    mse_rls(i) = 10*log10(compute_mse(x,xc_rls));
    mse_afa(i) = 10*log10(compute_mse(x,xc_afa));

    snr_lms(i) = compute_snr(x,xc_lms);
    snr_nlms(i) = compute_snr(x,xc_nlms);
    snr_rls(i) = compute_snr(x,xc_rls);
    snr_afa(i) = compute_snr(x,xc_afa);

    psnr_lms(i) = compute_psnr(x,xc_lms);
    psnr_nlms(i) = compute_psnr(x,xc_nlms);
    psnr_rls(i) = compute_psnr(x,xc_rls);
    psnr_afa(i) = compute_psnr(x,xc_afa);
end
% Compare SNR 
snr_before = compute_snr(x,xn);
snr_imp_lms = snr_lms - snr_before;
snr_imp_nlms = snr_nlms - snr_before;
snr_imp_rls = snr_rls - snr_before;
snr_imp_afa = snr_afa - snr_before;
figure;
plot(p,snr_imp_lms,p,snr_imp_nlms,p,snr_imp_rls,p,snr_imp_afa);
legend('LMS','NLMS','RLS','AFA');
xlabel('Filter Order');
ylabel('SNR_Imp (dB)');
title('Filter Order vs SNR Improvement');

figure;
plot(p,mse_lms,p,mse_nlms,p,mse_rls,p,mse_afa);
legend('LMS','NLMS','RLS','AFA');
xlabel('Filter Order');
ylabel('MSE');

figure;
plot(p,psnr_lms,p,psnr_nlms,p,psnr_rls,p,psnr_afa);
legend('LMS','NLMS','RLS','AFA');
xlabel('Filter Order');
ylabel('Peak SNR (dB)');





