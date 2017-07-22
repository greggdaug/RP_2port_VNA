clc
clf
clear all
close all

IP= '10.0.0.6';           % Input IP of your Red Pitaya...
port = 5000;
tcpipObj=tcpip(IP, port);
tcpipObj.InputBufferSize = 16384*32;
tcpipObj.OutputBufferSize = 16384*32;

%% Open connection with your Red Pitaya
fopen(tcpipObj);
tcpipObj.Terminator = 'CR/LF';
flushinput(tcpipObj)
flushoutput(tcpipObj)

%% Digital IO test
fprintf(tcpipObj,'DIG:PIN:DIR OUT,DIO7_P');
fprintf(tcpipObj,'DIG:PIN DIO7_P,0');
% fprintf(tcpipObj,'DIG:PIN? DIO7_P');
pause(0.1);
fprintf(tcpipObj,'DIG:PIN:DIR OUT,DIO6_P');
fprintf(tcpipObj,'DIG:PIN DIO6_P,0');
% fprintf(tcpipObj,'DIG:PIN? DIO6_P');
pause(0.1);  
fprintf(tcpipObj,'DIG:PIN:DIR OUT,DIO5_P');
fprintf(tcpipObj,'DIG:PIN DIO5_P,0');
% fprintf(tcpipObj,'DIG:PIN? DIO5_P');
pause(0.1);
fprintf(tcpipObj,'DIG:PIN:DIR OUT,DIO4_P');
fprintf(tcpipObj,'DIG:PIN DIO4_P,0');
% fprintf(tcpipObj,'DIG:PIN? DIO4_P');
pause(0.1);
fprintf(tcpipObj,'DIG:PIN:DIR OUT,DIO3_P');
fprintf(tcpipObj,'DIG:PIN DIO3_P,0');
% fprintf(tcpipObj,'DIG:PIN? DIO3_P');

%% Close connection with Red Pitaya
fclose(tcpipObj);