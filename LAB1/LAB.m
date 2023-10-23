clc
close all
clear 
A0 = 1;
f0 = 100;
fs = 10*f0;
fd= 100*fs;
tau = 10;
n1= 0:1/fd:3/f0 -1/fd;
n2= 0:1/fd:2/((2/3)*f0) -1/fd;
n3= 0:1/fd:1/((1/3)*f0) - 1/fd;
m1 = A0*cos(2*pi*f0*n1);
m2 = A0*square(2*pi*(2/3)*f0*n2);
m3 = A0*sawtooth(2*pi*(1/3)*f0*n3,1/2);
figure(1);
hold on
plot(n1,m1);
plot(n2,m2);
plot(n3,m3);
title('THE THREE SIGNALS');
xlabel('t(sec)');
ylabel('amplitude');
t = 0:  1/fd : (3/f0) - (1/fd);
clk =( A0*square(2*pi*fs*t,tau) + A0)/(2*A0);
clk(1501:1510)=1;   %there is a problem in the clock so i had to do this
clk(1901:1910)=1;   %there is a problem in the clock so i had to do this
figure(2);
plot(t,clk);
title('Pulses Train');
axis([0 3/f0-1/fd 0 2]);
xlabel('t(sec)');
ylabel('amplitude');
spikes = [diff(clk) 0];
comp = max(spikes ,0);
comp(1) = 1;
figure;
plot(comp);
axis([0 1000 -1 2]);
var1 = 0;
var2 = 0;
var3 = 0;
for i = 1:length(m1)
    if(comp(i) == 1)
        var1 = m1(i);
        var2 = m2(i);
        var3 = m3(i);
    end
    spam1(i) = var1;
    spam2(i) = var2;
    spam3(i) = var3;
    
end
pam1 = spam1.* clk;
pam2 = spam2.* clk;
pam3 = spam3.* clk;
figure(3);
subplot(2,2,1);
plot(pam1);
title('PAM FOR THE FIRST SIGNAL');
subplot(2,2,2);
plot(pam2);
title('PAM FOR THE SECOND SIGNAL');
subplot(2,2,3);
plot(pam3);
title('PAM FOR THE THIRD SIGNAL');
%%
% ............................. MUX .............................
T = [pam1(1)*ones(1,tau) zeros(1,10) pam2(1)*ones(1,tau) zeros(1,10) pam3(1)*ones(1,tau) zeros(1,10) ];
for i=101:100:3000
    if(clk(i) == 1) 
        
        T = [T pam1(i)*ones(1,tau) zeros(1,10) pam2(i)*ones(1,tau) zeros(1,10) pam3(i)*ones(1,tau) zeros(1,10)];
        
    end
    
end
figure(4);
plot(T);
title('THE RESULTED SIGNAL AFTER MULTIPLEXING ');
%%
%.............................DEMUX .............................
pam1_ret = zeros(1,length(pam1));
pam2_ret = zeros(1,length(pam2));
pam3_ret = zeros(1,length(pam3));
k = 0:100:length(pam1)-100;
d1 = 1;
d2 = 1;
d3 = 1;
for j=1:60:length(T)
  
   pam1_ret(1 + k(d1) :k(d1) + 10) = T(j)*ones(1,10);
   d1 = d1+1;
   
end
for j=21:60:length(T)
  
   pam2_ret(1 + k(d2) :k(d2) + 10) = T(j)*ones(1,10);
   d2 = d2+1;
   
end
for j=41:60:length(T)
  
   pam3_ret(1 + k(d3) :k(d3) + 10) = T(j)*ones(1,10);
   d3 = d3+1;
   
end
clc
figure(5);
plot(pam1_ret);
title('PAM FOR THE FIRST SIGNAL AFTER DEMULTIPLEXING');
figure(6);
plot(pam2_ret);
title('PAM FOR THE ESCOND SIGNAL AFTER DEMULTIPLEXING');
figure(7);
plot(pam3_ret);
title('PAM FOR THE THIRD SIGNAL AFTER DEMULTIPLEXING');
%%
%............................. FILTERING .............................

[b1,a1] =butter(4 ,4*f0/(fd/2)); 
sPAM1 = filter(b1,a1,pam1_ret);
[b2,a2] =butter(4 ,4*f0/(fd/2)); 
sPAM2 = filter(b2,a2,pam2_ret);
[b3,a3] =butter(4 ,4*f0/(fd/2)); 
sPAM3 = filter(b3,a3,pam3_ret);
hold on
plot(sPAM1);
plot(m1);
figure;
hold on
plot(sPAM2);
plot(m2);
figure;
hold on
plot(sPAM3);
plot(m3);
%%
%............................. SPECTRUM .............................
NP = length(pam1);
P1 = fft(pam1)/ NP;
P2 = fft(pam2)/NP;
P3 = fft(pam3)/NP;
PAM1 = fftshift(P1);
PAM2 = fftshift(P2);
PAM3 = fftshift(P3);
f = -fd/2 : fd/NP : fd/2 - fd/NP;
figure(8);
plot(f,abs(PAM1));
title('SPECTRUM OF FIRST PAM SIGNAL');
xlabel('FREQUENCY (Hz)');
ylabel('MAGNITUDE');
figure(9);
plot(f,abs(PAM2));
title('SPECTRUM OF SECOND PAM SIGNAL');
xlabel('FREQUENCY (Hz)');
ylabel('MAGNITUDE');
figure(10);
plot(f,abs(PAM3));
title('SPECTRUM OF THIRD PAM SIGNAL');
xlabel('FREQUENCY (Hz)');
ylabel('MAGNITUDE');
