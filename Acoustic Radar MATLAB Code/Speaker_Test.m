% Test Output of Speaker and Trasmit Line of Jammer System

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

SpeedSoundWave_ms = 343;           %[m/s]  ->speed of wave
Fc_Hz = 10e3;                      % [Hz]
Listening_Time = 5;                %[s]  
Fs = 44.1e3;                       %[Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
Ts = 1/Fs;                         % Sampling period
N = Listening_Time * Fs;           % number of samples
t = 0:Ts:(N-1)*Ts;                 % Time Vector for Pulse
frequencyaxis1 = (-N/2:1:(N/2-1))*Fs/N; % Frequency-axis to Plot Frequency Domain of Pulse

%% Record Received Samples from the Microphone

RecLength_s = Listening_Time; 
recObj = audiorecorder(Fs,24,1);
recordblocking(recObj, RecLength_s);   % Records audio for a fixed number of seconds
RX_signal = getaudiodata(recObj).';   % Store recorded audio signal in double-precision array

% Plot Signal Received from Speaker
figure(1); 
axes('fontsize', 12);

plot(t, RX_signal); % plot received signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (Linear)', 'fontsize', 12);
title('Speaker Signal', 'fontsize', 12);

% Plot Frequencies Received from Speaker
figure(2)
plot(frequencyaxis1,20*log10(abs(fftshift(fft(RX_signal)))));
title('Frequency Plot of Signal Received from Speaker');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
