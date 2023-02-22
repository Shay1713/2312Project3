%% KEY
clf;
close all

% creating audio recording parameters
Fs = 44100; 
nBits = 16; 
nChannels = 2; 
ID = -1;       % default audio input device 
load handel.mat % plotting parameters
window = hamming(512);
N_overlap = 256;
N_fft = 1024;


%% CITATIONS
% downsample: https://www.mathworks.com/help/signal/ref/downsample.html
% y = downsample(x, n); % n is downsampling factor
% y = downsample(x, n, phase); % where phase 0 to n-1

% upsample: https://www.mathworks.com/help/signal/ref/upsample.html
% y = upsample(x,n); 
% y = upsample(x,n,phase); 

% highpass: https://www.mathworks.com/help/signal/ref/highpass.html
% y = highpass(x,fpass,fs); 

% lowpass: https://www.mathworks.com/help/signal/ref/lowpass.html
% y = lowpass(x,fpass,fs); 


%% Code Professor Gave
%target_F = 8000
%sampling_freq = 44100

%stopband_st = target_F/sampling_freq
%passband_end = (target_F-2000)/sampling_freq

%F = [0 passband_end stopband_st 1] 
%A = [1 1 0 0] 
%lpf = firls(255, F, A); 
%filtered = filter(lpf, A, combined);


%% QUESTION 1 - design lp filter & downsampling implementation to only extracts speech content out

filename1 = 'phrase1.wav'; % phrase 1 from project 1
[x, Fs] = audioread(filename1); % calling the .wav file

% original phrase 1 recording
graph1 = graphing('Spectogram phrase1.wav', 'Time (seconds)', 'Frequency (Hz)', ...
     x, hamming(512), 256, 1024, 44100); % plotting original signal


% speech only phrase1.wav, then downsample and lowpass
% disp(size(x));
singleX = x(:,1); % reshape to be a single column
yDS = downsample(singleX, 2); % downsample
yDSLP = lowpass(yDS,3000,Fs); % lowpass PLAY AROUND WITH DIFFERENT CUTOFFS
%graph2 = graphing('Spectogram phrase1DSLP.wav', 'Time (seconds)', 'Frequency (Hz)', ...
    %yDSLP, hamming(512), 256, 1024, 44100); 
% sound(yDSLP);


%% QUESTION 2: downsample & lowpass & highpass filter

yLPFilterPreDS = lowpass(yDSLP,1000,Fs); % applying lowpass
% graph3 = graphing('Spectogram yLPFilterPreDS', 'Time (seconds)', 'Frequency (Hz)', ...
%     yLPFilterPreDS, hamming(512), 256, 1024, 44100);
yHPFilterPreDS = highpass(yDSLP,4000,Fs); % applying highpass
% graph4 = graphing('Spectogram yHPFilterPreDS', 'Time (seconds)', 'Frequency (Hz)', ...
%     yHPFilterPreDS, hamming(512), 256, 1024, 44100);

yLPFilterPostDS = downsample(yLPFilterPreDS, 2); % downsampling
% graph5 = graphing('Spectogram yLPFilterPostDS', 'Time (seconds)', 'Frequency (Hz)', ...
%     yLPFilterPostDS, hamming(512), 256, 1024, 44100);
yHPFilterPostDS = downsample(yHPFilterPreDS, 2); % downsampling
% graph6 = graphing('Spectogram yHPFilterPostDS', 'Time (seconds)', 'Frequency (Hz)', ...
%     yHPFilterPostDS, hamming(512), 256, 1024, 44100);



%% QUESTION 3: 

yLPFilterPreDS2 = lowpass(yLPFilterPostDS,1000,Fs); % applying lowpass
% graph7 = graphing('Spectogram yLPFilterPreDS2', 'Time (seconds)', 'Frequency (Hz)', ...
%      yLPFilterPreDS2, hamming(512), 256, 1024, 44100);
yHPFilterPreDS2 = highpass(yHPFilterPostDS,4000,Fs); % applying highpass
% graph8 = graphing('Spectogram yHPFilterPreDS2', 'Time (seconds)', 'Frequency (Hz)', ...
%     yHPFilterPreDS2, hamming(512), 256, 1024, 44100);

yLPFilterPostDS2 = downsample(yLPFilterPreDS2, 2); % apply downsample
% graph9 = graphing('Spectogram yLPFilterPostDS2', 'Time (seconds)', 'Frequency (Hz)', ...
%     yLPFilterPostDS2, hamming(512), 256, 1024, 44100);
yHPFilterPostDS2 = downsample(yHPFilterPreDS2, 2); % apply lowpass
% graph10 = graphing('Spectogram yHPFilterPostDS2', 'Time (seconds)', 'Frequency (Hz)', ...
%      yHPFilterPostDS2, hamming(512), 256, 1024, 44100);




%% QUESTION 4

% reconstruct the signal by applying everything but opposite
yLPFilterPostDS21 = upsample(yLPFilterPostDS2, 2); % upsample reverse 1a
yLPFilterPostDS211 = lowpass(yLPFilterPostDS21,1000,Fs); % lowpass reverse 1
yHPFilterPostDS21 = upsample(yHPFilterPostDS2, 2); % upsample reverse 1b
yHPFilterPostDS211 = highpass(yHPFilterPostDS21,4000,Fs); % highpass reverse1
 
yLPFilterPostDS1 = upsample(yLPFilterPostDS211, 2); % upsample reverse 2a
yLPFilterPostDS11 = lowpass(yLPFilterPostDS1,1000,Fs); % lowpass reverse 2
yHPFilterPostDS1 = upsample(yHPFilterPostDS211, 2); % upsample reverse 2b
yHPFilterPostDS11 = highpass(yHPFilterPostDS1,4000,Fs); % upsample reverse 2
synthSig11 = [yLPFilterPostDS1 , yHPFilterPostDS1]; % this is recombining the signals

synthSig1 = upsample(synthSig11, 2); % reversing intial downsample
synthSig = lowpass(synthSig1,1000,Fs); % matching intial lowpass

% save as 'team[8]-synthesized.wav'
[synthesis, Fs] = writeReadFile('team[8]-synthesized.wav', synthSig, Fs); % saving to a .wav file

% playing the resulting sound and describe it
% sound(synthesis); % what sound like

% plot spectrogram
% graph11 = graphing('Spectogram team[8]-synthesized.wav', 'Time (seconds)', 'Frequency (Hz)', ...
%      synthesis, hamming(512), 256, 1024, 44100); % plotting synthesized file




%% FUNCTION STORAGE ********************

function [x, Fs] = writeReadFile(fileNameDo, x, Fs) % READ AND WRITE AUDIO FILES 
    filename1 = fileNameDo; % creating a .wav file for the sine tone
    audiowrite(filename1,x,Fs); % writing and reading the sine tone into the .wav file
    [x,Fs] = audioread(filename1);
end


function fig = graphing(plotTitle, xTitle, yTitle, file, window, N_overlap, N_fft, Fs) % PLOTTING SPECTROGRAMS
    [S, F, T, P] = spectrogram(file(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); 
    fig = figure;
    surf(T, F, 10*log10(P), 'edgecolor', 'none');
    axis tight;
    view(0,90);
    colormap(jet);
    set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
    xlabel(xTitle); ylabel(yTitle); title(plotTitle); % labels
end

