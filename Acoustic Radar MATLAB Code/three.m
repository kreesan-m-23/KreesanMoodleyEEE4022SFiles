% Audio Radar Project - Lowpass Filtering After Down-mixing

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

SpeedSoundWave_ms = 343;           %[m/s]  ->speed of wave
Fc = 10e3;                      % [Hz]
PulseWidth_s = 1;                  %[s]  
Fs = 44.1e3;                       %[Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
Ts = 1/Fs;                         % Sampling period
N = PulseWidth_s * Fs;             % number of samples
t = 0:Ts:(N-1)*Ts;                 % Time Vector for Pulse
frequencyaxis1 = (-N/2:1:((N/2)-1))*(Fs/N); % Frequency-axis to Plot Frequency Domain of Pulse

%% Generate the transmit pulse

TxPulse = sin(2*pi*(8000*t + 2000*t.^2));

%% Exercise 3 - downmix the transmit pulse

% Multiply by cos(x) and -sin(x) to get I and Q channels, respectively

I = TxPulse.*cos(2*pi*Fc*t);
Q = TxPulse.*-sin(2*pi*Fc*t);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);
    % LPF with Fpass = 2100 Hz, Fstop = 6000 Hz, Apass = .001 dB, Astop = 80 dB
    % Made for specific bandwidth of 4 kHz

FilterCoeff = [-1.53785105157198e-05,0.000130395417810422,0.000307205296087191,0.000493641064532186,0.000523071237085720,0.000206486459690951,-0.000559895414633393,-0.00166139944430918,-0.00268699258188829,-0.00298940878748828,-0.00192092207848543,0.000801932573902799,0.00474654559300812,0.00862180669894579,0.0105096233804997,0.00851018316305392,0.00165013391372272,-0.00927040324639787,-0.0212870916561839,-0.0296423227547979,-0.0289728325141533,-0.0150430984819859,0.0135234716553144,0.0541621450120427,0.100431503965337,0.143358564544868,0.173741781831962,0.184715552228408,0.173741781831962,0.143358564544868,0.100431503965337,0.0541621450120427,0.0135234716553144,-0.0150430984819859,-0.0289728325141533,-0.0296423227547979,-0.0212870916561839,-0.00927040324639787,0.00165013391372272,0.00851018316305392,0.0105096233804997,0.00862180669894579,0.00474654559300812,0.000801932573902799,-0.00192092207848543,-0.00298940878748828,-0.00268699258188829,-0.00166139944430918,-0.000559895414633393,0.000206486459690951,0.000523071237085720,0.000493641064532186,0.000307205296087191,0.000130395417810422,-1.53785105157198e-05];
 
    Signal_I_AferLPF = filter(FilterCoeff,1,I);
    Signal_Q_AferLPF = filter(FilterCoeff,1,Q);

output = complex(Signal_I_AferLPF, Signal_Q_AferLPF);

% Plot the down-mixed signal in the frequency domain

figure(1);
axes('fontsize', 12);
plot(frequencyaxis1,(abs(fftshift(fft(output)))));
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Amplitude', 'fontsize', 12);
title('Output: Frequency-domain', 'fontsize', 12);
