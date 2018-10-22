% Audio Radar Project - Generate Chirp Signal, Trasmit and Record Signals and Bandpass Filter

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

SpeedSoundWave_ms = 343;           %[m/s]  ->speed of wave
Fc_Hz = 10e3;                      % [Hz]
PulseWidth_s = 3;                  %[s]  
Fs = 44.1e3;                       %[Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
Ts = 1/Fs;                         % Sampling period
N = PulseWidth_s * Fs;             % number of samples
t = 0:Ts:(N-1)*Ts;                 % Time Vector for Pulse
frequencyaxis1 = (-N/2:1:(N/2-1))*Fs/N; % Frequency-axis to Plot Frequency Domain of Pulse

%% Generate the transmit pulse

% Pure sinusoid 
TxPulse = sin(2*pi*10000*t + 2*pi*2000*t.^2);
Tx_Signal_Single = sin(2*pi*(10000*(t - 1/2) + 1*4000*(t - 1/2).^2) ).* rect( (t - 1/2)/1 )   ;  % transmit signal with chirp pulse

%% Create the transmit signal

NumOfZeros = N;
TransmitSignal = [zeros(1, NumOfZeros) Tx_Signal_Single zeros(1, NumOfZeros)]; % -> Create the transmit signal

TimeAxis_TxSignal_s = (0:1:(length(TransmitSignal)-1))*Ts;

%% Play out transmit signal through the speakers

soundsc(TransmitSignal,Fs, 24) % Transmit the signal

%% Record received samples from the microphone

RecLength_samples = length(TransmitSignal);
RecLength_s = RecLength_samples*1/Fs; 
recObj = audiorecorder(Fs,24,1);
recordblocking(recObj, RecLength_s);  % Records audio for a fixed number of seconds
RX_signal1 = getaudiodata(recObj).';   % Store recorded audio signal in double-precision array

% BPF with Fs = 44100 Hz, Fstop1 = 5500, Fpass1 = 7900, Fpass2 = 12100, Fstop2 = 14500, Astop1 = Astop2 = 40 dB, Apass = .001 dB 
FilterCoeff = [-0.000652609550284060,-0.00226080014835753,0.000282954093856288,-0.00291998509081353,0.000207102506938406,0.00415100128354289,0.000624789500247805,-0.000900631818761759,0.00148118374194383,-0.00555423246018557,-0.00726682717122072,0.00657763557982094,0.00724617553716979,-0.000539184462538904,0.00747653097660437,-0.000916743341238398,-0.0242223124735308,-0.00621822583281655,0.0174812607683140,0.00226944324231889,0.0126464654678576,0.0279905466974503,-0.0307911717290630,-0.0539587747623631,0.0144610502385423,0.0156455657754238,0.000487803111235825,0.100426767847059,0.0471070810924698,-0.219462352600266,-0.157155233688628,0.243124041091914,0.243124041091914,-0.157155233688628,-0.219462352600266,0.0471070810924698,0.100426767847059,0.000487803111235825,0.0156455657754238,0.0144610502385423,-0.0539587747623631,-0.0307911717290630,0.0279905466974503,0.0126464654678576,0.00226944324231889,0.0174812607683140,-0.00621822583281655,-0.0242223124735308,-0.000916743341238398,0.00747653097660437,-0.000539184462538904,0.00724617553716979,0.00657763557982094,-0.00726682717122072,-0.00555423246018557,0.00148118374194383,-0.000900631818761759,0.000624789500247805,0.00415100128354289,0.000207102506938406,-0.00291998509081353,0.000282954093856288,-0.00226080014835753,-0.000652609550284060];
RX_signal = filter(FilterCoeff,1,RX_signal1);

figure(1); 
axes('fontsize', 12);

subplot(2,1,1);
plot(TimeAxis_TxSignal_s, TransmitSignal); % plot received signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit signal', 'fontsize', 12);
grid on;

subplot(2,1,2);
plot(TimeAxis_TxSignal_s, RX_signal); % plot received signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Received signal', 'fontsize', 12);
grid on;

figure(2)
frequencyaxis1 = (-length(RX_signal)/2:1:(length(RX_signal)/2-1))*Fs/length(RX_signal);
subplot(2,1,1);
plot(frequencyaxis1,(abs(fftshift(fft(TransmitSignal)))));
title('Frequency Plot of Transmit Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude (Linear)');
subplot(2,1,2);
plot(frequencyaxis1,20*log10(abs(fftshift(fft(RX_signal)))));
title('Frequency Plot of Received Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
