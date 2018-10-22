% Audio Radar Project - Down-mixing of Chirp Signal

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

SpeedSoundWave_ms = 343;           %[m/s]  ->speed of wave
Fc = 10e3;                         % [Hz]
PulseWidth_s = 1;                  %[s]  
Fs = 44.1e3;                       %[Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
Ts = 1/Fs;                         % Sampling period
N = PulseWidth_s * Fs;             % number of samples
t = 0:Ts:(N-1)*Ts;                 % Time Vector for Pulse
frequencyaxis1 = (-N/2:1:((N/2)-1))*(Fs/N); % Frequency-axis to Plot Frequency Domain of Pulse

%% Generate the transmit pulse

TxPulse = sin(2*pi*(8000*t + 2000*t.^2));

%% Exercise 2 - downmix the transmit pulse 

% Multiply by cos(2*pi*fc*t) and -sin(2*pi*fc*t) to get I and Q channels, respectively

I = TxPulse.*cos(2*pi*Fc*t);
Q = TxPulse.*-sin(2*pi*Fc*t);

% Generate the output, which are complex numbers. The I and Q parts are real
% part. Tip: output = I + 1i*Q; % or complex(I, Q)

output = complex(I,Q);

% Plot the transmit and down-mixed signals in the frequency domain
figure(1);
axes('fontsize', 12);

subplot(2,1,1)
plot(frequencyaxis1,(abs(fftshift(fft(TxPulse)))));
xlabel('Frequency (Hz', 'fontsize', 12);
ylabel('Amplitude (Linear)', 'fontsize', 12);
title('Original Pulse: Frequency-domain', 'fontsize', 12);

subplot(2,1,2)
plot(frequencyaxis1,(abs(fftshift(fft(output)))));
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Amplitude (Linear)', 'fontsize', 12);
title('Down-mixed Pulse: Frequency-domain', 'fontsize', 12);
