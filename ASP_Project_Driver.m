%% ASP Term Project Driver
% Adaptive Noise Cancelation in Music signals
% Chris Novak, Natalie Meyer

%% Read in audio data
audiodir = './ASP_Project_Audio/';
listname = dir(audiodir);
listname = listname(3:end);
fs = 44100;
t_per_song = 10; % 20 second clips of each song
num_samples = t_per_song * fs;
music_files = {};
for i = 1:length(listname)
    [y,fs] = audioread([audiodir, listname(i).name],[1 num_samples]);
    music_files{i} = y;
end


%% Add Noise

x = music_files{i};
% make mono for now
x = mean(x,2);
xlen = length(x);
gain = .01;
mean = 0;
noise = gain * randn(xlen,1) + mean;
figure;

subplot(311);
plot(x);
xnoise = x + noise;
subplot(312);
plot(noise); 
subplot(313);
plot(xnoise);

%% Pre process WOLA




