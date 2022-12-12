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
 t_per_song = 5; % 10 second clips of each song
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

[xn,ref_noise] = create_and_add_noise(x,.5,10,.5,'crowd');
% xlen = length(x);
% noiselen = 1 * fs; % change noise every 1 sec
% num_frames = xlen / noiselen;
% xn = zeros(xlen,1); % signal + noise
% ref_noise = zeros(xlen,1); % noise
% gains = [.05 .05 .05 .05 .05 .05 .05 .05 .05 .05]; % can change this to be varying
% for f = 0:num_frames-1
%     ind = 1+f*noiselen:noiselen+f*noiselen;
%     ref_noise(ind) = gains(f+1)*randn(noiselen,1);
% end
% 
% % form some correlation between primary and reference noise
% lp = fir1(10,.4);
% prim_noise = filter(lp,1,ref_noise); 
% xn = x + prim_noise;

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
sigma = 1; % initial update matrix param
gamma = .5; % gain parameter for afa convergence between .5 and 1

xc_lms = perform_lms(xn,ref_noise,mu,p);
xc_nlms = perform_nlms(xn,ref_noise,mu,p);
xc_rls = perform_rls(xn,ref_noise,lambda,sigma,p);
xc_afa = perform_afa(xn,ref_noise,gamma,p);

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
ylabel('Mean Square Error');
title('LMS mu = 0.001')
xlabel('Samples n (iterations)')
subplot(412);
plot(converge_nlms);
ylabel('Mean Square Error');
title('NLMS mu = 0.001')
xlabel('Samples n (iterations)')
subplot(413);
plot(converge_rls);
subplot(414);
plot(converge_afa);

Compare SNR 
snr_before = compute_snr(x,xn)

snr_lms = compute_snr(x,xc_lms)

snr_nlms = compute_snr(x,xc_nlms)

snr_rls = compute_snr(x,xc_rls)

snr_afa = compute_snr(x,xc_afa)

% Compare PSNR 
psnr_before = compute_psnr(x,xn)

psnr_lms = compute_psnr(x,xc_lms)

psnr_nlms = compute_psnr(x,xc_nlms)

psnr_rls = compute_psnr(x,xc_rls)

psnr_afa = compute_psnr(x,xc_afa)

%% Multirate

xc_wavelet = perform_wavelet_anc(xn,ref_noise,.01,10,1,'db6');

snr_mr = compute_snr(x,xc_wavelet)

%% Test parameters
order = 4:20;
for p = 1:length(order)
    mu = linspace(.0001,.1,100);
    mse_lms = zeros(1,length(mu));
    mse_nlms = zeros(1,length(mu));
    for i = 1:length(mu)
        xc_lms = perform_lms(xn,ref_noise,mu(i),order(p));
        xc_nlms = perform_nlms(xn,ref_noise,mu(i),order(p));
    
        mse_lms(i) = compute_mse(x,xc_lms);
        mse_nlms(i) = compute_mse(x,xc_nlms);
    end
%     
%     figure;
%     plot(mu,mse_lms,mu,mse_nlms);
%     xlabel('Step Size (mu)');
%     ylabel('MSE');
%     legend('LMS','NLMS');
    
    [min_mse_lms,loc_lms] = min(mse_lms);
    best_mu_lms(p) = mu(loc_lms);
    [min_mse_nlms,loc_nlms] = min(mse_nlms);
    best_mu_nlms(p) = mu(loc_nlms);
    
    lambda = linspace(.98,1,10);
    mse_rls = zeros(1,length(lambda));
    for i = 1:length(lambda)
        xc_rls = perform_rls(xn,ref_noise,lambda(i),1,order(p));
    
        mse_rls(i) = compute_mse(x,xc_rls);
    end
%     figure;
%     plot(lambda,mse_rls);
%     xlabel('Forgetting Factor (lambda)');
%     ylabel('MSE');
%     legend('RLS');
     [min_mse_rls,loc_rls] = min(mse_rls);
     best_lam_rls(p) = lambda(loc_rls);
    
    gamma = linspace(.5,1,50);
    mse_afa = zeros(1,length(gamma));
    for i = 1:length(gamma)
        xc_afa = perform_afa(xn,ref_noise,gamma(i),order(p));
    
        mse_afa(i) = compute_mse(x,xc_afa);
    end
%     figure;
%     plot(gamma,mse_afa);
%     xlabel('Gain Step Size (gamma)');
%     ylabel('MSE');
%     legend('AFA');
    [min_mse_afa,loc_afa] = min(mse_afa);
    best_gam_afa(p) = gamma(loc_afa);
end

figure;
plot(order,best_mu_lms(1:17),order,best_mu_nlms(1:17));
xlabel('Filter Order');
ylabel('Mean Square Error');
title('Best Mu vs Filter Order');
legend('LMS','NLMS');


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
    

%% 
converge_lms = abs(x - xc_lms);
converge_nlms = abs(x - xc_nlms);
converge_rls = abs( x - xc_rls);
converge_afa = abs(x - xc_afa);
figure;
sgtitle('Convergence of Algorithms with Optimal Parameters')

subplot(411);
plot(converge_lms);
title('LMS')
subplot(412);
plot(converge_nlms);
title('NMLS')
subplot(413);
plot(converge_rls);
title('RLS')
subplot(414);
plot(converge_afa);
title('AFA')

snr_before = compute_snr(x,xn)

snr_lms = compute_snr(x,xc_lms)

snr_nlms = compute_snr(x,xc_nlms)

snr_rls = compute_snr(x,xc_rls)

snr_afa = compute_snr(x,xc_afa)





