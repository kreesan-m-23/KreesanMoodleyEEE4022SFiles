% Test Output of Microphone of Jammer System

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

SpeedSoundWave_ms = 343;           %[m/s]  ->speed of wave
Fc_Hz = 9.5e3;                     % [Hz]
Trasmit_Time = 3;                  %[s]  
Fs = 44.1e3;                       %[Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
Ts = 1/Fs;                         % Sampling period
N = Trasmit_Time * Fs;             % number of samples
t = 0:Ts:(N-1)*Ts;                 % Time Vector for Pulse
frequencyaxis1 = (-N/2:1:(N/2-1))*Fs/N; % Frequency-axis to Plot Frequency Domain of Pulse
%% Generate the transmit pulse

% Pure sinusoid 
TxPulse = sin(2*pi*10000*t + 2*pi*2000*t.^2);
Tx_Signal_Single = sin(2*pi*(Fc_Hz*(t - 1/2) + 0*0*(t - 1/2).^2) ).* rect( (t - Trasmit_Time/2)/Trasmit_Time )   ;  % transmit signal

%% Create the transmit signal

NumOfZeros = N;
TransmitSignal = Tx_Signal_Single;% -> Create the transmit signal

TimeAxis_TxSignal_s = (0:1:(length(TransmitSignal)-1))*Ts;

%% Play out transmit signal through the speakers

soundsc(TransmitSignal,Fs, 16) % Transmit the signal
