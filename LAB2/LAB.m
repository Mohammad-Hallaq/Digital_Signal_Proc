clear;clc;close all;
Fs=1e4;% Sampling Frequency
bps=16;% Number of bits per sample
Dur = 3; %recording duration
record = audiorecorder(Fs,bps,1);% Audio Recorder Object
disp('Please Start Recording:');% Display a message to user
recordblocking(record, Dur);% Record Sound for 3 Seconds
disp('End of Recording.')
play(record);
x = getaudiodata(record);% Get sound data from recObject
Ts = 1/Fs;% Sampling Period
T = length(x);% Get number of total samples from recObject
n = 0 : Ts : Dur - Ts;% Define Disctere Time axis
plot(n,x),grid on;% Plot the signal
xlabel('Time (sec)'); ylabel('Amplitude');% Set x, y axis names
title('Speech Signal');% Set a title for the plot
audiowrite('soundtest.wav',x,Fs);% Saving the record
%%
[out_rec , Fs] = audioread('verygoodsound.wav');% Reading the outer record
disp(Fs);
figure;
plot(n,out_rec),grid on;% Plot the signal
xlabel('Time (sec)'); ylabel('Amplitude');% Set x, y axis names
title('Speech Signal');% Set a title for the plot
A_let = out_rec(9165 : 14800);% Getting A letter
sound(A_let);
figure;
plot(n(1 : 14800 - 9165 +1),A_let),grid on;% Plot the signal
xlabel('Time (sec)'); ylabel('Amplitude');% Set x, y axis names
title('A letter');% Set a title for the plot
audiowrite('A.wav',A_let,Fs);% Saving the A voice
win_A = A_let(3000:3199);
figure;
plot(n(1:1199 - 1000 +1),win_A),grid on;% Plot the signal
xlabel('Time (sec)'); ylabel('Amplitude');% Set x, y axis names
title('A window of letter A ');% Set a title for the plot
N = 256;
A_ftr = fft(win_A,N);
A_spectrum = fftshift(A_ftr);
f = -Fs/2 : Fs/N : Fs/2 - Fs/N;
figure;
plot(f,abs(A_spectrum));grid on;% Plot the signal
xlabel('Frequency (Hz)'); ylabel('Amplitude');% Set x, y axis names
title('Fourier Transform for "A" letter window ');% Set a title for the plot
%%
P = 12;
[a,g] = lpc(win_A,P);
freq = 0 : Fs/N : Fs/2 - Fs/N;
padding_pars = [a zeros(1,length(win_A)-13)]';
PAR_FOUR = fft(padding_pars,N);
figure;
plot(freq ,-20*log10( abs(PAR_FOUR(1:128))),freq ,20*log10(abs(A_ftr(1:128))),'r'); 
xlabel('Frequency (Hz)'); ylabel('Amplitude (dB)');% Set x, y axis names
title('Comparison between Spectrum and Estimation of letter A window');
%%
[h,w] = freqz(1,a,128);
sig = var(win_A);
eps = sqrt(sig)*randn(1,N);%% The Noise ~ N(0,sig)
s = filter(12,a,eps);
[pxx,F] = pcov(s,12);
figure;
plot(freq,20*log10(abs(h)))
title('PSD estimatin for "A" letter');
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
hold on
plot(freq,10*log10(pxx(1:128)))
legend('LPC PSD Estimate','pcov PSD Estimate')
%%
Blocks_num = length(out_rec)/200;% Number of Energy Blocks 
Energy = zeros(1,Blocks_num);% Energy of Windows
j=0;
for i=1:150
    c=out_rec(200*j+1:200*(j+1));
    Energy(1,i)=sum(c.^2);
    j=j+1;
end
t_e = 0:3/Blocks_num:3-3/Blocks_num; 
figure;
plot(n,out_rec),grid on;% Plot the signal
xlabel('Time (sec)'); ylabel('Amplitude');% Set x, y axis names
hold on
plot(t_e,Energy/max(Energy),'r')

