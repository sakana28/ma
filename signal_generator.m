%% Simulated localized and distributed fault in rolling
% element bearing
%
% G. D’Elia and M. Cocconcelli
clear
clc
%% Bearing geometry
d = 21.4; % bearing roller diameter [mm]
D = 203; % pitch circle diameter [mm]
n = 23; % number of rolling elements
contactAngle = 9*pi/180; % contact angle
faultType = 'inner';

%% Speed profile
N = 2048; % number of points per revolution
Ltheta = 10000*N; % signal length
theta = (0:Ltheta-1)*2*pi/N;
fc = 10;
fd = 0.08*fc;
fm = 0.1*fc;
fr = fc + 2*pi*fd.*(cumsum(cos(fm.*theta)/N));

%% Localized fault
varianceFactor = 0.04;
fs = 20000; % sample frequency [Hz]
k = 2e13;
zita = 5/100;
fn = 6e3; % natural frequency [Hz]
Lsdof = 2^8;
SNR_dB = 0;
qAmpMod = 0.3;
 [tLocal,xLocal,xNoiseLocal,frTimeLocal,meanDeltaTLocal,varDeltaTLocal,meanDeltaTimpOverLocal, ... 
     varDeltaTimpOverLocal,errorDeltaTimpLocal] = bearingSignalModelLocal(d,D,contactAngle,n, ...
     faultType,fr,fc,fd,fm,N,varianceFactor,fs,k,zita,fn,Lsdof,SNR_dB,qAmpMod);

%% Distributed fault
fs = 20000; % sample frequency [Hz]
SNR_dB = 0;
qFault = 1;
qStiffness = 0.1;
qRotation = 0.1;
[tDist,xDist,xNoiseDist,frTimeDist] = bearingSignalModelDist(d,D,contactAngle,n,faultType,fc,fd,fm ...
    ,fr,N,fs,SNR_dB,qFault,qStiffness,qRotation);

%% Draw distributed fault

figure(1)
plot(tDist,xDist);
xlabel('time [s]');
ylabel('amplitude');
xlim([0 0.5]);
ylim([-11 11]);
title('Distributed fault (inner) without noise')

figure(2)
plot(tDist,xNoiseDist);
xlabel('time [s]');
ylabel('amplitude');
xlim([0 0.5]);
ylim([-11 11]);
title('Distributed fault (inner) with noise')


    
    