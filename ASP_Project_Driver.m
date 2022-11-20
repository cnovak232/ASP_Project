%% ASP Term Project Driver
% Adaptive Noise Cancelation in Music signals
% Chris Novak, Natalie Meyer

%% Read in audio data
audiodir = './ASP_Project_Audio/';
listname = dir(audiodir);
listname = listname(3:end);
fs = 44100;
t_per_song = 10; % 5 second clips of each song
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
sn = zeros(xlen,1); % noise
gains = ones(1,noiselen)*.05; % can change this 

for f = 0:num_frames-1
    ind = 1+f*noiselen:noiselen+f*noiselen;
    noise = gains(f+1)*randn(noiselen,1);
    xn(ind) = x(ind) + noise;
    sn(ind) = noise;
end

figure;
subplot(311);
plot(x);
subplot(312);
plot(sn); 
subplot(313);
plot(xn);

%% Run adaptive noise cancelation process

% win_len = 512; % define in samples for now
% 
% % hann window - 50% overlap negates need for normalization
% win = hann(win_len,"periodic")';
% hop = win_len / 2; 
% 
% ys = run_anc(xn,sn,win,hop);

xc = perform_lms(xn,sn,.1,20);

figure;
subplot(311);
plot(x);
subplot(312);
plot(xn);
subplot(313);
plot(xc);

%% Compare SNR








