% Audio Radar Project - Matched Filtering in Frequency-domain

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
N = length(t);          % Number of Samples
RangeLineAxis_m = t*c/2;% Range-axis
frequencyaxis1 = (-N/2:1:((N/2)-1))*(fs/N); % Frequency-axis to Plot Frequency Domain of Pulse

% Target parameters 
R_target = 300;         % Target range 
td = 2* R_target/c;     % Two way time-delay 

R_target2 = 500;        % Target range 
td2 = 2* R_target2/c;   % Two way time-delay 

%% Create the transmit signal

Tx_Signal = cos(2*pi*(fc*(t - T/2) + 0.5*K*(t - T/2).^2) ).* rect( (t - T/2)/T )   ;  % transmit signal with chirp pulse

%% Generate the transmit pulse
 
NumSamplesTxPulse = ceil(T/ts);             % number of samples of the transmit pulse 
Tx_p = Tx_Signal(1: NumSamplesTxPulse);     % transmit pulse only

%% Generate a simulated received signal - not transmitting as yet

% Generate the Received Signal
 
Rx_Signal = (cos(2*pi*(fc*(t - T/2 - td + T/2) + 0.5*K*(t - T/2 - td + T/2).^2) ).* rect( (t - T/2 - td + T/2)/T )) + (cos(2*pi*(fc*(t - T/2 - td2 + T/2) + 0.5*K*(t - T/2 - td2 + T/2).^2) ).* rect( (t - T/2 - td2 + T/2)/T ));  % received signal   
% extra T/2 to get peak of target at correct range

% Plot Transmit and Simulated Received Signals
figure(1)
axes('fontsize', 12);

subplot(2,1,1);
plot(t, Tx_Signal);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (Linear)', 'fontsize', 12);
title('Transmit Pulse', 'fontsize', 12);
grid on;
 
RangeLineAxis_m = t*c/2;

subplot(2,1,2);
plot( RangeLineAxis_m,Rx_Signal);
xlabel('Range (m)', 'fontsize', 12);
ylabel('Amplitude (Linear)', 'fontsize', 12); 
title('Simulated Range Line', 'fontsize', 12);
grid on;

%% Exercise 2-3: downmixing

I_Tx = Tx_Signal.*cos(2*pi*fc*t);
Q_Tx = Tx_Signal.*-sin(2*pi*fc*t);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % LPF with Fpass = 2100 Hz, Fstop = 6000 Hz, Apass = .001 dB, Astop = 80 dB
    % Made for specific bandwidth of 4 kHz

FilterCoeff = [-1.53785105157198e-05,0.000130395417810422,0.000307205296087191,0.000493641064532186,0.000523071237085720,0.000206486459690951,-0.000559895414633393,-0.00166139944430918,-0.00268699258188829,-0.00298940878748828,-0.00192092207848543,0.000801932573902799,0.00474654559300812,0.00862180669894579,0.0105096233804997,0.00851018316305392,0.00165013391372272,-0.00927040324639787,-0.0212870916561839,-0.0296423227547979,-0.0289728325141533,-0.0150430984819859,0.0135234716553144,0.0541621450120427,0.100431503965337,0.143358564544868,0.173741781831962,0.184715552228408,0.173741781831962,0.143358564544868,0.100431503965337,0.0541621450120427,0.0135234716553144,-0.0150430984819859,-0.0289728325141533,-0.0296423227547979,-0.0212870916561839,-0.00927040324639787,0.00165013391372272,0.00851018316305392,0.0105096233804997,0.00862180669894579,0.00474654559300812,0.000801932573902799,-0.00192092207848543,-0.00298940878748828,-0.00268699258188829,-0.00166139944430918,-0.000559895414633393,0.000206486459690951,0.000523071237085720,0.000493641064532186,0.000307205296087191,0.000130395417810422,-1.53785105157198e-05];

    Signal_I_Tx_AferLPF = filter(FilterCoeff,1,I_Tx);
    Signal_Q_Tx_AferLPF = filter(FilterCoeff,1,Q_Tx);

% Generate complex numbers with I being real part and Q being the imaginary
% part. Tip: x = complex(a,b);

Tx_Baseband = complex(Signal_I_Tx_AferLPF, Signal_Q_Tx_AferLPF);

I_Rx = Rx_Signal.*cos(2*pi*fc*t);
Q_Rx = Rx_Signal.*-sin(2*pi*fc*t);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);
   
    Signal_I_Rx_AferLPF = filter(FilterCoeff,1,I_Rx);
    Signal_Q_Rx_AferLPF = filter(FilterCoeff,1,Q_Rx);

% Generate complex numbers with I being real part and Q being the imaginary
% part. Tip: x = complex(a,b);

Rx_Baseband = complex(Signal_I_Rx_AferLPF, Signal_Q_Rx_AferLPF);

figure(2);
axes('fontsize', 12);
plot(frequencyaxis1,((abs(fftshift(fft(Tx_Baseband))))));
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Amplitude', 'fontsize', 12);
title('Baseband Output Pulse: Frequency-domain', 'fontsize', 12);

%% Exercise 4: Pulse compression or matched filtering

% Obtain the matched filter

Matched_Filter = conj(fliplr(Tx_Baseband));

% Obtain pulse compression output: convolution between the matched filter
% and the baseband received signal

output = ifft(fft(Rx_Baseband) .* fft(Matched_Filter));

% Plot the pulse compressed output  

figure(3)

t_add = 0:ts:(T/2); % To account for time delay and include 0-0.5 sec in time range
out = [zeros(1,length(t_add)) output];
o_new = abs(out(1:length(RangeLineAxis_m))); % fixing length of vector
plot(RangeLineAxis_m,o_new);
title('Matched Filter Output Showing Peaks at Ranges');
xlabel ('Range (m)');
ylabel('Amplitude(Linear)');
grid On;
