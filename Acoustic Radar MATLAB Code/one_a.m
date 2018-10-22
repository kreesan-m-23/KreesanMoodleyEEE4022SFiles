% Audio Radar Project - Define Parameters and Plot Transmit Pulse

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

SpeedSoundWave_ms = 343;           % [m/s] Speed of Wave
Fc = 10e3;                         % [Hz] Centre Frequency of Pulse
PulseWidth_s = 1;                  % [s]  Length of Plot
Fs = 44.1e3;                       % [Hz] Sampling Rate
Ts = 1/Fs;                         % [s] Sampling Period
N = PulseWidth_s * Fs;             % number of samples
t = 0:Ts:(N-1)*Ts;                 % Time Vector for Pulse
frequencyaxis1 = (-N/2:1:(N/2-1))*Fs/N; % Frequency-axis to Plot Frequency Domain of Pulse

%% Generate the transmit pulse

TxPulse = cos(2*pi*Fc*t); % Sinusoidal Pulse

%% Plot the transmit signal

% Time-domain Plot of Signal
figure(1);
axes('fontsize', 12);
plot(t,TxPulse);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit Pulse: Time-domain', 'fontsize', 12);

% Frequency-domain Plot of Signal
figure(2);
axes('fontsize', 12);
plot(frequencyaxis1,(abs((fft(TxPulse)))));
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Amplitude', 'fontsize', 12);
title('Transmit Pulse: Frequency-domain', 'fontsize', 12);
