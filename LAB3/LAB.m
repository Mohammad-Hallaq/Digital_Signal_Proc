clc
close all
clear 
Fs = 16000;
nbits = 8;
channels = 1;
DUR = 5;
%%
recorder = audiorecorder(Fs,nbits,channels);
disp('Start speaking.')
recordblocking(recorder, DUR);
disp('End of Recording.');
play(recorder);
%%
[out_rec , Fs] = audioread('snd.wav');
t = 0:1/Fs : 1 - 1/Fs;
%myRecording = getaudiodata(out_rec);
figure;
plot(out_rec);
figure;
plot(t,out_rec);
S_let = out_rec(1300:2688);
O_let = out_rec(2751:4017);
figure;
plot(S_let);
figure;
plot(O_let);
sound(S_let);
sound(O_let);
%%
FOU_s = fft(S_let);
FOU_o = fft(O_let);
FOU_S = fftshift(FOU_s);
FOU_O = fftshift(FOU_o);
Ns = length(S_let);
fs = -Fs/2 : Fs/Ns : Fs/2 - Fs/Ns;
No = length(O_let);
fo = -Fs/2 : Fs/No : Fs/2 - Fs/No;
figure;
plot(fs,abs(FOU_S));
figure;
plot(fo,abs(FOU_O));
%%
Ns = 512;
No = 512;
FOU_s = fft(S_let,Ns);
FOU_o = fft(O_let,No);
FOU_S = fftshift(FOU_s);
FOU_O = fftshift(FOU_o);
fs = -Fs/2 : Fs/Ns : Fs/2 - Fs/Ns;
fo = -Fs/2 : Fs/No : Fs/2 - Fs/No;
figure;
plot(fs,abs(FOU_S));
figure;
plot(fo,abs(FOU_O));
