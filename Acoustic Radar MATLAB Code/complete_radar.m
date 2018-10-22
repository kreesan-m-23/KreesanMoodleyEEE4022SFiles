% Audio Radar Project - Working Acoustic Radar Programme

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

c = 343;                           % Speed of sound wave [m/s]

% Radar parameters
fc = 10e3;                         % [Hz]
T = 0.03;                          % Pulse length [s] 
fs = 44.1e3;                       % Sampling rate in [Hz]. So 44 100 samples are obtained in one second
ts = 1/fs;                         % Sampling period [s]
B = 4e3;                           % Bandwidth [Hz]
PRI = .06;                         % Pulse repitition interval [s]
PRF = 1/PRI;                       % [Hz]
NumPulses = 32;                    % Number of pulses 
R_max = c*PRI/2;                   % Maximum range of target to simulate [m]
t_max = PRI*NumPulses;             % Maximum time to simulate 
K = B/T;                           % Chirp rate
lamda = c/fc;                      % Wavelength of transmitted pulse
PulseCompressionGain = T*B;        % Pulse gain due to bandwidth and pulse length
BlindRange_m = c*T/2/PulseCompressionGain;   % Range which is too close to be detected
UnambiguousRange_m = c*PRI/2;                % Maximum range that can be unambiguously detected [m]
RangeResolution_m = c/(2*B);                 % Minimum separation needed between targets to detect both individually
UnambiguousVelocity_cms = PRF/2*lamda/2*100; % Maximum velocity that can be unambiguously detected [m]

% Vectors and axes for calculations and plotting
t  = (0:ts:(t_max-ts));                      % Time vector (in seconds) 
N_total = length(t);                         % Number of samples in full-time vector
N_single_line = length(t)/NumPulses;         % Number of samples in one PRI
time_single_line = t(1):ts:t(N_single_line); % Time vector for one PRI
RangeLineAxis_m = time_single_line*c/2;      % Range line axis
frequencyaxis1 = (-N_single_line/2:1:((N_single_line/2)-1))*(fs/N_single_line); % Frequency axis for single PRI
frequencyaxis2 = (-N_total/2:1:((N_total/2)-1))*(fs/N_total);                   % Frequency axis to plot full-time signal
DopplerHz = (-NumPulses/2:1:(NumPulses/2-1))*PRF/NumPulses; % Doppler frequency axis
Velocity_cmpersec = DopplerHz*(lamda/2)*100;                % Velocity axis derived from Doppler frequency [cm/s]

% Accommodate for 2 Extra Pulses Sent 
Num_Pulses_Extra = 2;
time_extra = Num_Pulses_Extra * PRI;
extra_vector_length = time_extra * fs; % Number of samples in extra vector
extra_vector = zeros(1, extra_vector_length);
tmax_with_extra_vector = t_max + time_extra;
time_total_with_extra = (0:ts:(tmax_with_extra_vector));

%% Display radar parameters

clc
disp(' ');
disp(['Blind Range = ' num2str(roundn(BlindRange_m,-3)) ' m']);
disp(['Unambiguous Range = ' num2str(roundn(UnambiguousRange_m,0)) ' m']);
disp(['Range Resolution = ' num2str(roundn(RangeResolution_m,-3)) ' m']);
disp(['Unambiguous Velocity = ' num2str(roundn(UnambiguousVelocity_cms,-2)) ' cm/s']);
disp(['Pulse Compression gain = ' num2str(roundn(PulseCompressionGain,0)) '']);
disp(' ');

%% Generate the Transmit Signal

Tx_Signal_Single = sin(2*pi*(fc*(time_single_line - T/2) + 0.5*K*(time_single_line - T/2).^2)).* rect((time_single_line - T/2)/T );  % transmit signal with chirp pulse
Tx_Signal = repmat(Tx_Signal_Single,1,NumPulses);
Tx_Signal = [extra_vector Tx_Signal];
x = Tx_Signal_Single';
soundsc(Tx_Signal,fs, 24) % Transmit the signal

%% Capture Echoes

recObj = audiorecorder(fs,24,1);
recordblocking(recObj, (tmax_with_extra_vector + time_extra));  % Records audio for a fixed number of seconds
Rx_Signal1 = getaudiodata(recObj).';   % Store recorded audio signal in double-precision array

Tx_Signal = Tx_Signal((length(extra_vector)+1):length(Tx_Signal));
Rx_Signal1 = Rx_Signal1((length(extra_vector)+1):(length(Rx_Signal1)- length(extra_vector)));

% BPF with Fs = 44100 Hz, Fstop1 = 5500, Fpass1 = 7900, Fpass2 = 12100, Fstop2 = 14500, Astop1 = Astop2 = 40 dB, Apass = .001 dB 
FilterCoeff_BPF = [-0.000652609550284060,-0.00226080014835753,0.000282954093856288,-0.00291998509081353,0.000207102506938406,0.00415100128354289,0.000624789500247805,-0.000900631818761759,0.00148118374194383,-0.00555423246018557,-0.00726682717122072,0.00657763557982094,0.00724617553716979,-0.000539184462538904,0.00747653097660437,-0.000916743341238398,-0.0242223124735308,-0.00621822583281655,0.0174812607683140,0.00226944324231889,0.0126464654678576,0.0279905466974503,-0.0307911717290630,-0.0539587747623631,0.0144610502385423,0.0156455657754238,0.000487803111235825,0.100426767847059,0.0471070810924698,-0.219462352600266,-0.157155233688628,0.243124041091914,0.243124041091914,-0.157155233688628,-0.219462352600266,0.0471070810924698,0.100426767847059,0.000487803111235825,0.0156455657754238,0.0144610502385423,-0.0539587747623631,-0.0307911717290630,0.0279905466974503,0.0126464654678576,0.00226944324231889,0.0174812607683140,-0.00621822583281655,-0.0242223124735308,-0.000916743341238398,0.00747653097660437,-0.000539184462538904,0.00724617553716979,0.00657763557982094,-0.00726682717122072,-0.00555423246018557,0.00148118374194383,-0.000900631818761759,0.000624789500247805,0.00415100128354289,0.000207102506938406,-0.00291998509081353,0.000282954093856288,-0.00226080014835753,-0.000652609550284060];
Rx_Signal = filter(FilterCoeff_BPF,1,Rx_Signal1);

% Generate matrix of range-lines
RxSignalMatrix1 = reshape(Rx_Signal,  length(Rx_Signal)/NumPulses, NumPulses);
RxSignalMatrix = RxSignalMatrix1.'; % Each row is a range line

%% Down-mixing

% LPF with Fpass = 2100 Hz, Fstop = 6000 Hz, Apass = .001 dB, Astop = 80 dB
% Made for specific bandwidth of 4 kHz
FilterCoeff_LPF = [-1.53785105157198e-05,0.000130395417810422,0.000307205296087191,0.000493641064532186,0.000523071237085720,0.000206486459690951,-0.000559895414633393,-0.00166139944430918,-0.00268699258188829,-0.00298940878748828,-0.00192092207848543,0.000801932573902799,0.00474654559300812,0.00862180669894579,0.0105096233804997,0.00851018316305392,0.00165013391372272,-0.00927040324639787,-0.0212870916561839,-0.0296423227547979,-0.0289728325141533,-0.0150430984819859,0.0135234716553144,0.0541621450120427,0.100431503965337,0.143358564544868,0.173741781831962,0.184715552228408,0.173741781831962,0.143358564544868,0.100431503965337,0.0541621450120427,0.0135234716553144,-0.0150430984819859,-0.0289728325141533,-0.0296423227547979,-0.0212870916561839,-0.00927040324639787,0.00165013391372272,0.00851018316305392,0.0105096233804997,0.00862180669894579,0.00474654559300812,0.000801932573902799,-0.00192092207848543,-0.00298940878748828,-0.00268699258188829,-0.00166139944430918,-0.000559895414633393,0.000206486459690951,0.000523071237085720,0.000493641064532186,0.000307205296087191,0.000130395417810422,-1.53785105157198e-05];

I_Tx = Tx_Signal_Single.*cos(2*pi*fc*time_single_line);
Q_Tx = Tx_Signal_Single.*-sin(2*pi*fc*time_single_line);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);

    Signal_I_Tx_AferLPF = filter(FilterCoeff_LPF,1,I_Tx);
    Signal_Q_Tx_AferLPF = filter(FilterCoeff_LPF,1,Q_Tx);

% Generate complex numbers with I real and Q being complex

Tx_Baseband = complex(Signal_I_Tx_AferLPF, Signal_Q_Tx_AferLPF);

I_Rx = Rx_Signal.*cos(2*pi*fc*t);
Q_Rx = Rx_Signal.*-sin(2*pi*fc*t);

% Low pass filter the I and Q channel.
    % Use fdatool to calculate the coefficients of the LPF 
    % Example code:  Signal_I_AferLPF = filter(FilterCoeff,1,InputSignal_I);
    % Example code:  Signal_Q_AferLPF = filter(FilterCoeff,1,InputSignal_Q);

    Signal_I_Rx_AferLPF = filter(FilterCoeff_LPF,1,I_Rx);
    Signal_Q_Rx_AferLPF = filter(FilterCoeff_LPF,1,Q_Rx);

% Generate complex numbers with I real and Q being complex

Rx_Baseband = complex(Signal_I_Rx_AferLPF, Signal_Q_Rx_AferLPF);
Rx_Baseband_Matrix = (reshape(Rx_Baseband,  length(Rx_Baseband)/NumPulses, NumPulses)).';

%% Matched Filtering

Tx_Baseband_Modified = [Tx_Baseband zeros(1,length(Rx_Baseband)-length(Tx_Baseband))];
H = fft(conj(fliplr(Tx_Baseband_Modified)));
MF1 = ifft(fft(Rx_Baseband,[],2).*H,[],2);

MF = reshape(MF1,  length(MF1)/NumPulses, NumPulses).';

%% Generate window matrix

Window = hamming(NumPulses);
numcols = size(Rx_Baseband_Matrix,2);
Window_Matrix = repmat(Window,1,numcols);

%% Generate Range-Doppler Map

% Generate Range-Doppler matrix by Doppler-processing
output1 = abs((fftshift(fft(Window_Matrix.*MF,[],1),1)));

% Left-Circular Shift Output to 'Zero' Signals Coming Directly from Speaker
index_of_max = find(output1 == max(output1(:)));  % Find index of maximum value
col_of_max = ceil(index_of_max/NumPulses);        % Find column in which highest value is
output = circshift(output1,[0,-(col_of_max -1)]); % Circular shift left

% Left-Circular Shift the Range Lines to 'Zero' Signals Coming Directly from Speaker
index_of_max = find(MF == max(MF(:)));     % Find index of maximum value
col_of_max = ceil(index_of_max/NumPulses); % Find column in which highest value is
MF = circshift(MF,[0,-(col_of_max -1)]);   % Circular shift left

%% Plot Signal Results

figure(1)
subplot(2,1,1);
plot(t,Tx_Signal);
xlabel('Time (s)');
ylabel('Magnitude');
title('Time Plot of Transmit Signal');

subplot(2,1,2);
plot(t,abs(MF1));
xlabel('Time');
ylabel('Magnitude');
title('Received Matched-Filtered Signal');

figure(2)
subplot(2,1,1);
plot(frequencyaxis2,20*log10(abs(fftshift(fft(Rx_Baseband)))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Frequency Spectrum of Baseband Received Signal');

subplot(2,1,2);
plot(frequencyaxis2,20*log10(abs(fftshift(fft(MF1)))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Frequency Spectrum of Matched Filtering Output');

figure(3)
plot(RangeLineAxis_m, abs(MF(1,:)));
xlabel('Range (m)');
ylabel('Amplitude');
title('Range Line of First Received Pulse');

figure(4)
subplot(2,1,1);
plot(frequencyaxis2,20*log(abs(fftshift(fft(Tx_Signal)))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Frequency Spectrum of Transmit Signal');

subplot(2,1,2);
plot(frequencyaxis2,abs(fftshift(fft(Rx_Signal))));
xlabel('Frequency');
ylabel('Magnitude (dB)');
title('Frequency Spectrum of Received Signal');

figure(5)
imagesc(RangeLineAxis_m, 1:NumPulses, abs(MF));
xlabel('Range (m)');
ylabel('Pulse Number');
title('Plot of Rangle-lines Matrix');

% Plot Range-Doppler map
figure(6)
imagesc(RangeLineAxis_m(:,200:length(output)-1000), Velocity_cmpersec, output(:,200:length(output)-1000));
colorbar;
grid On;
ax = gca;
ax.GridColor = [1 1 1];
colormap('jet');
axis xy;
xlabel('Range (m)');
ylabel('Velocity (cm/s)');
title('Range-Doppler Map');
