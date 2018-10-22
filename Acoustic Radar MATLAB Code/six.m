% Audio Radar Project - Generate Matrix Plot of Range Lines

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

c = 343;                % [m/s]  ->speed of wave
 
% Radar parameters
fc = 10e3;              % [Hz]
T = .25;                % Pulse length in [s]  
fs = 44.1e3;            % [Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
ts = 1/fs;              % Sampling period
B = 4e3;                % Bandwidth (Hz) 
NumPulses = 8;          % Number of pulses 
PRI = 1;                % Pulse Repetition Interval, Ru = 600
PRF = 1/PRI;
 
t_max = PRI;            % Maximum time to simulate 
K = B/T;                % Chirp rate
t  = (0:ts:(t_max-ts)); % Time vector (in seconds) 

time = []; % Creation of complete time vector
for i = 0:NumPulses-1
    time = [time (t+i*t(length(t)))];
end
 
% Target parameters 
R_target = 50;                       % Target range 
Vel_target = 0.003;                  % Target radial velocity in m/s
td = 2 * R_target/c;                 % Time taken for signal to travel to target and back to the radar

N = length(t);
frequencyaxis1 = (-N/2:1:((N/2)-1))*(fs/N);
RangeLineAxis_m = t*c/2;

%% Generate the transmit pulse

Tx_Signal_Single = cos(2*pi*(fc*(t - T/2) + 0.5*K*(t - T/2).^2) ).* rect( (t - T/2)/T )   ;  % transmit signal with chirp pulse

%% Create the transmit signal

Tx_Signal = repmat(Tx_Signal_Single,1,NumPulses);

%% Generate a simulated received signal - not transmitting as yet
%{
Rx_Signal_Single = (cos(2*pi*(fc*(t - T/2 - td + T/2) + 0.5*K*(t - T/2 - td + T/2).^2) ).* rect( (t - T/2 - td + T/2)/T )) + (cos(2*pi*(fc*(t - T/2 - td2 + T/2) + 0.5*K*(t - T/2 - td2 + T/2).^2) ).* rect( (t - T/2 - td2 + T/2)/T ));  % received signal   
Rx_Signal = repmat(Rx_Signal_Single,1,NumPulses);
%}
Rx_Signal1 = zeros(1, size(time,2));
 
 for Count_PulseNum = 1: NumPulses
 
    tdn = 2*(R_target - Vel_target*PRI*(Count_PulseNum - 1))/c + PRI*(Count_PulseNum - 1); 
     
    Rx_Signal1 = Rx_Signal1 +  cos(2*pi*(fc*(time - T/2 - tdn) + 0.5*K*(time - T/2 - tdn).^2) ).* rect( (time - T/2 - tdn)/T );  % received signal   
 
 end
 
% BPF with Fs = 44100 Hz, Fstop1 = 5500, Fpass1 = 7900, Fpass2 = 12100, Fstop2 = 14500, Astop1 = Astop2 = 40 dB, Apass = .001 dB 
FilterCoeff_BPF = [-0.000652609550284060,-0.00226080014835753,0.000282954093856288,-0.00291998509081353,0.000207102506938406,0.00415100128354289,0.000624789500247805,-0.000900631818761759,0.00148118374194383,-0.00555423246018557,-0.00726682717122072,0.00657763557982094,0.00724617553716979,-0.000539184462538904,0.00747653097660437,-0.000916743341238398,-0.0242223124735308,-0.00621822583281655,0.0174812607683140,0.00226944324231889,0.0126464654678576,0.0279905466974503,-0.0307911717290630,-0.0539587747623631,0.0144610502385423,0.0156455657754238,0.000487803111235825,0.100426767847059,0.0471070810924698,-0.219462352600266,-0.157155233688628,0.243124041091914,0.243124041091914,-0.157155233688628,-0.219462352600266,0.0471070810924698,0.100426767847059,0.000487803111235825,0.0156455657754238,0.0144610502385423,-0.0539587747623631,-0.0307911717290630,0.0279905466974503,0.0126464654678576,0.00226944324231889,0.0174812607683140,-0.00621822583281655,-0.0242223124735308,-0.000916743341238398,0.00747653097660437,-0.000539184462538904,0.00724617553716979,0.00657763557982094,-0.00726682717122072,-0.00555423246018557,0.00148118374194383,-0.000900631818761759,0.000624789500247805,0.00415100128354289,0.000207102506938406,-0.00291998509081353,0.000282954093856288,-0.00226080014835753,-0.000652609550284060];
Rx_Signal = filter(FilterCoeff_BPF,1,Rx_Signal1);
% This will have to be adjusted to methods shown in 1.b) when transmitting

%% Exercise 6: Generate the range line for multiple pulses, assuming a stationary target
    % Tips: use reshape to convert a vector into a matrix
    % use imagesc() to plot a matrix. 

Rx_Signal_Matrix1 = reshape(Rx_Signal,  length(Rx_Signal)/NumPulses, NumPulses);
Rx_Signal_Matrix = Rx_Signal_Matrix1.'; % Each row is a range line

%% Exercise 2-3: downmixing
% LPF with Fpass = 2100 Hz, Fstop = 6000 Hz, Apass = .001 dB, Astop = 80 dB
% Made for specific bandwidth of 4 kHz
FilterCoeff_LPF = [-1.53785105157198e-05,0.000130395417810422,0.000307205296087191,0.000493641064532186,0.000523071237085720,0.000206486459690951,-0.000559895414633393,-0.00166139944430918,-0.00268699258188829,-0.00298940878748828,-0.00192092207848543,0.000801932573902799,0.00474654559300812,0.00862180669894579,0.0105096233804997,0.00851018316305392,0.00165013391372272,-0.00927040324639787,-0.0212870916561839,-0.0296423227547979,-0.0289728325141533,-0.0150430984819859,0.0135234716553144,0.0541621450120427,0.100431503965337,0.143358564544868,0.173741781831962,0.184715552228408,0.173741781831962,0.143358564544868,0.100431503965337,0.0541621450120427,0.0135234716553144,-0.0150430984819859,-0.0289728325141533,-0.0296423227547979,-0.0212870916561839,-0.00927040324639787,0.00165013391372272,0.00851018316305392,0.0105096233804997,0.00862180669894579,0.00474654559300812,0.000801932573902799,-0.00192092207848543,-0.00298940878748828,-0.00268699258188829,-0.00166139944430918,-0.000559895414633393,0.000206486459690951,0.000523071237085720,0.000493641064532186,0.000307205296087191,0.000130395417810422,-1.53785105157198e-05];

I_Tx = Tx_Signal_Single.*cos(2*pi*fc*t);
Q_Tx = Tx_Signal_Single.*-sin(2*pi*fc*t);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);

    Signal_I_Tx_AferLPF = filter(FilterCoeff_LPF,1,I_Tx);
    Signal_Q_Tx_AferLPF = filter(FilterCoeff_LPF,1,Q_Tx);

% Generate complex numbers with I being real part and Q being the imaginary
% part. Tip: x = complex(a,b);

Tx_Baseband = complex(Signal_I_Tx_AferLPF, Signal_Q_Tx_AferLPF);

I_Rx = Rx_Signal.*cos(2*pi*fc*time);
Q_Rx = Rx_Signal.*-sin(2*pi*fc*time);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);

    Signal_I_Rx_AferLPF = filter(FilterCoeff_LPF,1,I_Rx);
    Signal_Q_Rx_AferLPF = filter(FilterCoeff_LPF,1,Q_Rx);

% Generate complex numbers with I being real part and Q being the imaginary
% part. Tip: x = complex(a,b);

Rx_Baseband = complex(Signal_I_Rx_AferLPF, Signal_Q_Rx_AferLPF);

Rx_Baseband_Matrix1 = reshape(Rx_Baseband,  length(Rx_Baseband)/NumPulses, NumPulses);
Rx_Baseband_Matrix = Rx_Baseband_Matrix1.'; % Each row is a range line

%% Exercise 4/6: Pulse compression or matched filtering

Tx_Baseband_Modified = [Tx_Baseband zeros(1,length(Rx_Baseband)-length(Tx_Baseband))];
H = fft(conj(fliplr(Tx_Baseband_Modified)));
MF1 = ifft(fft(Rx_Baseband,[],2).*H,[],2);
MF = reshape(MF1,  length(MF1)/NumPulses, NumPulses).'; 

% Plot Range Lines in Matrix
figure(1)
imagesc(RangeLineAxis_m, (1:1:NumPulses), abs(MF));
colorbar;
grid On;
ax = gca;
ax.GridColor = [1 1 1];
colormap('jet');
axis xy;
xlabel('Range (m)');
ylabel('Pulse Number');
title('Range Map');

% Plot Specific Range Lines 1 and 'NumPulses' to Show How they Become Bright on the
% Range Map
figure(2)

subplot(2,1,1)
plot(RangeLineAxis_m, abs(MF(1,:)));
title('First Range Line of Matrix');
xlabel('Range (m)');
ylabel('Amplitude (Linear)');

subplot(2,1,2)
plot(RangeLineAxis_m, abs(MF(NumPulses,:)));
title('Last Range Line of Matrix');
xlabel('Range (m)');
ylabel('Amplitude (Linear)');
