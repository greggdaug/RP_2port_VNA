clc
clf
clear all
close all

fstart=3e6;
fstop=10e6;
N=16383;
fs=125e6;
ts=1/fs;
T=ts*N;
fbin=fs/N;
freq=0:fs/N:fs/2;
index_fstart=round(fstart/fbin);
index_fstop=round(fstop/fbin);
k=(fstop-fstart)/T;
t=0:ts:N*ts;
x1=sin(2*pi*(fstart.*t+k./2*t.^2));

% IP= '10.0.0.6';           % Input IP of your Red Pitaya...
% port = 5000;
% tcpipObj=tcpip(IP, port);
% tcpipObj.InputBufferSize = 16384*32;
% tcpipObj.OutputBufferSize = 16384*32;

vna_digio_set(1,1,0);

%% Read & plot
[signal_num1,signal_num2]=rp_getdata(x1);
c1=fft(signal_num1);
c2=fft(signal_num2);
z1=c2./-c1;
figure(1);
plot(signal_num1)
figure(2);
plot(signal_num2)
figure(3);
plot(freq(index_fstart:index_fstop)/1e6,20*log10(abs(z1(index_fstart:index_fstop))),'b')
% hold on
grid on

%% Close connection with Red Pitaya
fclose(tcpipObj);