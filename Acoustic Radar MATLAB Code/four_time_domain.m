% Audio Radar Project - Matched Filtering in Time-domain

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

c = 343;                % [m/s]  ->speed of wave
 
% Radar parameters
fc = 10e3;              % [Hz]
T = .03;                % Pulse length in [s]  
fs = 44.1e3;            % [Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
ts = 1/fs;              % Sampling period
B = 4e3;                % Bandwidth (Hz) 
R_max = 20;             % Maximum range of target to simulate
 
t_max = 2*R_max/c;      % Maximum time to simulate 
K = B/T;                % Chirp rate
t  = (0:ts:(t_max-ts)); % Time vector (in seconds) 
N = length(t);          % Number of Samples
RangeLineAxis_m = t*c/2;% Range-axis
frequencyaxis1 = (-N/2:1:((N/2)-1))*(fs/N); % Frequency-axis to Plot Frequency Domain of Pulse

% Target parameters 
R_target = 3;         % Target range 
td = 2* R_target/c;   % Two way time-delay 

R_target2 = 10;       % Target range 
td2 = 2* R_target2/c; % Two way time-delay 

%% Create the transmit signal

Tx_Signal = cos(2*pi*(fc*(t - T/2) + 0.5*K*(t - T/2).^2) ).* rect( (t - T/2)/T )   ;  % transmit signal with chirp pulse

%% Generate a simulated received signal - not transmitting as yet

Rx_Signal = (cos(2*pi*(fc*(t - T/2 - td + T/2) + 0.5*K*(t - T/2 - td + T/2).^2) ).* rect( (t - T/2 - td + T/2)/T )) + (cos(2*pi*(fc*(t - T/2 - td2 + T/2) + 0.5*K*(t - T/2 - td2 + T/2).^2) ).* rect( (t - T/2 - td2 + T/2)/T ));  % received signal      
% extra T/2 to get peak of target at correct range

% BPF with Fs = 44100 Hz, Fstop1 = 5500, Fpass1 = 7900, Fpass2 = 12100, Fstop2 = 14500, Astop1 = Astop2 = 40 dB, Apass = .001 dB 
FilterCoeff_BPF = [-0.000652609550284060,-0.00226080014835753,0.000282954093856288,-0.00291998509081353,0.000207102506938406,0.00415100128354289,0.000624789500247805,-0.000900631818761759,0.00148118374194383,-0.00555423246018557,-0.00726682717122072,0.00657763557982094,0.00724617553716979,-0.000539184462538904,0.00747653097660437,-0.000916743341238398,-0.0242223124735308,-0.00621822583281655,0.0174812607683140,0.00226944324231889,0.0126464654678576,0.0279905466974503,-0.0307911717290630,-0.0539587747623631,0.0144610502385423,0.0156455657754238,0.000487803111235825,0.100426767847059,0.0471070810924698,-0.219462352600266,-0.157155233688628,0.243124041091914,0.243124041091914,-0.157155233688628,-0.219462352600266,0.0471070810924698,0.100426767847059,0.000487803111235825,0.0156455657754238,0.0144610502385423,-0.0539587747623631,-0.0307911717290630,0.0279905466974503,0.0126464654678576,0.00226944324231889,0.0174812607683140,-0.00621822583281655,-0.0242223124735308,-0.000916743341238398,0.00747653097660437,-0.000539184462538904,0.00724617553716979,0.00657763557982094,-0.00726682717122072,-0.00555423246018557,0.00148118374194383,-0.000900631818761759,0.000624789500247805,0.00415100128354289,0.000207102506938406,-0.00291998509081353,0.000282954093856288,-0.00226080014835753,-0.000652609550284060];
Rx_Signal = filter(FilterCoeff_BPF,1,Rx_Signal);

% Plot Transmit and Simulated Received Signals
figure(1)
axes('fontsize', 12);

subplot(2,1,1);
plot(t, Tx_Signal);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit Pulse', 'fontsize', 12);
grid on;

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
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);
    % LPF with Fpass = 2100 Hz, Fstop = 6000 Hz, Apass = .001 dB, Astop = 80 dB
    % Made for specific bandwidth of 4 kHz

FilterCoeff_LPF = [-1.53785105157198e-05,0.000130395417810422,0.000307205296087191,0.000493641064532186,0.000523071237085720,0.000206486459690951,-0.000559895414633393,-0.00166139944430918,-0.00268699258188829,-0.00298940878748828,-0.00192092207848543,0.000801932573902799,0.00474654559300812,0.00862180669894579,0.0105096233804997,0.00851018316305392,0.00165013391372272,-0.00927040324639787,-0.0212870916561839,-0.0296423227547979,-0.0289728325141533,-0.0150430984819859,0.0135234716553144,0.0541621450120427,0.100431503965337,0.143358564544868,0.173741781831962,0.184715552228408,0.173741781831962,0.143358564544868,0.100431503965337,0.0541621450120427,0.0135234716553144,-0.0150430984819859,-0.0289728325141533,-0.0296423227547979,-0.0212870916561839,-0.00927040324639787,0.00165013391372272,0.00851018316305392,0.0105096233804997,0.00862180669894579,0.00474654559300812,0.000801932573902799,-0.00192092207848543,-0.00298940878748828,-0.00268699258188829,-0.00166139944430918,-0.000559895414633393,0.000206486459690951,0.000523071237085720,0.000493641064532186,0.000307205296087191,0.000130395417810422,-1.53785105157198e-05];

    Signal_I_Tx_AferLPF = filter(FilterCoeff_LPF,1,I_Tx);
    Signal_Q_Tx_AferLPF = filter(FilterCoeff_LPF,1,Q_Tx);

% Generate complex numbers with I being real part and Q being the imaginary
% part. Tip: x = complex(a,b);

Tx_Baseband = complex(Signal_I_Tx_AferLPF, Signal_Q_Tx_AferLPF);

I_Rx = Rx_Signal.*cos(2*pi*fc*t);
Q_Rx = Rx_Signal.*-sin(2*pi*fc*t);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);

    Signal_I_Rx_AferLPF = filter(FilterCoeff_LPF,1,I_Rx);
    Signal_Q_Rx_AferLPF = filter(FilterCoeff_LPF,1,Q_Rx);

% Generate complex numbers with I being real part and Q being the imaginary
% part. Tip: x = complex(a,b);

Rx_Baseband = complex(Signal_I_Rx_AferLPF, Signal_Q_Rx_AferLPF);

figure(2);
axes('fontsize', 12);
plot(frequencyaxis1,((abs(fftshift(fft(Tx_Baseband))))));
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Amplitude', 'fontsize', 12);
title('Output Pulse: Frequency-Domain', 'fontsize', 12);

%% Exercise 4: Pulse compression or matched filtering

% Obtain the matched filter

Matched_Filter = conj(fliplr(Tx_Baseband));

% Obtain pulse compression output: convolution between the matched filter
% and the baseband received signal

Output = conv(Matched_Filter, Rx_Baseband);

% Plot the pulse compressed output  

figure(3)

t_add = 0:ts:T/2; % To account for time delay and include 0-0.5 sec in time range
out = Output((length(t)-length(t_add)+1):length(Output)-1); % real vaules of matched filtering output and eliminating unimportant convolution values
o_new = out(1:length(RangeLineAxis_m)); % fixing length of vector

plot(RangeLineAxis_m,abs(o_new));
xlabel('Range (m)');
ylabel('Amplitude (Linear)');
title('Matched Filter Output Showing Peaks at Ranges');
grid On;
