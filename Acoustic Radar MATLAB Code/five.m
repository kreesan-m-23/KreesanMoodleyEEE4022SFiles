% Audio Radar Project - Create Multiple Pulses

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

c = 343;                % [m/s]  ->speed of wave
 
% Radar parameters
fc = 10e3;              % [Hz]
T = 1;                  % Pulse length in [s]  
fs = 44.1e3;            % [Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
ts = 1/fs;              % Sampling period
B = 4e3;                % Bandwidth (Hz) 
R_max = 600;            % Maximum range of target to simulate
 
t_max = 2*R_max/c;      % Maximum time to simulate 
K = B/T;                % Chirp rate
t  = (0:ts:(t_max-ts)); % Time vector (in seconds) 
 
% Target parameters 
R_target = 300;         % Target range 
td = 2* R_target/c;     % Two way time-delay 

N = length(t);
frequencyaxis1 = (-N/2:1:((N/2)-1))*(fs/N);

%% Generate the transmit pulse

Tx_Signal_Single = cos(2*pi*(fc*(t - T/2) + 0.5*K*(t - T/2).^2) ).* rect( (t - T/2)/T )   ;  % transmit signal with chirp pulse

%% Create the transmit signal for N = 8 pulses 

P = 8; % Number of Pulses

Tx_Signal = repmat(Tx_Signal_Single,1,P);

time = [];
for i = 0:P-1
    time = [time (t+i*t(length(t)))];
end

%% Plot transmit and simulated received signal

figure(1)
axes('fontsize', 12);

subplot(2,1,1)
plot(time, [Tx_Signal_Single zeros(1,(length(Tx_Signal)-length(Tx_Signal_Single)))]);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (Linear)', 'fontsize', 12);
title('Transmit Signal with Single Pulse', 'fontsize', 12);
ylim([-2 2])
grid on;

subplot(2,1,2)
plot(time, Tx_Signal);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (Linear)', 'fontsize', 12);
title('Transmit Signal with Multiple Pulses', 'fontsize', 12);
ylim([-2 2])
grid on;
