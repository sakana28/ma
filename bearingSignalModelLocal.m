function [t,x,xNoise,frTime,meanDeltaT,varDeltaT,meanDeltaTimpOver,varDeltaTimpOver,errorDeltaTimp] = ...
    bearingSignalModelLocal(d,D,contactAngle,n,faultType,fr,fc,fd,fm,N,varianceFactor,fs,k,zita,fn,Lsdof,SNR_dB,qAmpMod)
%% Generation of a simulated signal for localized fault in rolling element bearing
%
% Input:
% d = bearing roller diameter [mm]
% D = pitch circle diameter [mm]
% contactAngle = contact angle [rad]
% n = number of rolling elements
% faultType = fault type selection: inner, outer, ball [string]
% fr = row vector containing the rotation frequency profile
% fc = row vector containing the carrier component of the speed
% fm = row vector containing the modulation frequency
% fd = row vector containing the frequency deviation
% N = number of points per revolution
% varianceFactor = variance for the generation of the random contribution (ex. 0.04)
% fs = sample frequency of the time vector
% k = SDOF spring stiffness [N/m]
% zita = SDOF damping coefficient
% fn = SDOF natural frequency [Hz]
% Lsdof = length of the in number of points of the SDOF response
% SNR_dB = signal to noise ratio [dB]
% qAmpMod = amplitude modulation due to the load (ex. 0.3)
%
% Output:
% t = time signal [s]
% x = simulated bearing signal without noise
% xNoise = simulated bearing signal with noise
% frTime = speed profile in the time domain [Hz]
% meanDeltaT = theoretical mean of the inter-arrival times
% varDeltaT = theoretical variance of the inter-arrival times
% menDeltaTimpOver = real mean of the inter-arrival times
% varDeltaTimpOver = real variance of the inter-arrival times
% errorDeltaTimp = generated error in the inter-arrival times
%
% G. D’Elia and M. Cocconcelli

    if nargin < 14
        qAmpMod = 1;
    end

    switch faultType
        case 'inner'
            geometryParameter = 1 / 2 * (1 + d/D*cos(contactAngle)); % inner race fault
        case 'outer'
            geometryParameter = 1 / 2 * (1 - d/D*cos(contactAngle)); % outer race fault
        case 'ball'
            geometryParameter = 1 / (2*n) * (1 - (d/D*cos(contactAngle))^2)/(d/D); % outer race fault
    end

    Ltheta = length(fr);
    theta = (0:Ltheta-1)*2*pi/N;

    deltaThetaFault = 2*pi/(n*geometryParameter);
    numberOfImpulses = floor(theta(end)/deltaThetaFault);
    meanDeltaTheta = deltaThetaFault;
    varDeltaTheta = (varianceFactor*meanDeltaTheta)^2;
    deltaThetaFault = sqrt(varDeltaTheta)*randn([1 numberOfImpulses-1]) + meanDeltaTheta;
    thetaFault = [0 cumsum(deltaThetaFault)];
    frThetaFault = interp1(theta,fr,thetaFault,'spline');
    deltaTimp = deltaThetaFault ./ (2*pi*frThetaFault(2:end));
    tTimp = [0 cumsum(deltaTimp)];

    L = floor(tTimp(end)*fs); % signal length
    t = (0:L-1)/fs;
    frTime = interp1(tTimp,frThetaFault,t,'spline');

    deltaTimpIndex = round(deltaTimp*fs);
    errorDeltaTimp = deltaTimpIndex/fs - deltaTimp;

    indexImpulses = [1 cumsum(deltaTimpIndex)];
    index = length(indexImpulses);
    while indexImpulses(index)/fs > t(end)
        index = index - 1;
    end
    indexImpulses = indexImpulses(1:index);

    meanDeltaT = mean(deltaTimp);
    varDeltaT = var(deltaTimp);
    meanDeltaTimpOver = mean(deltaTimpIndex/fs);
    varDeltaTimpOver = var(deltaTimpIndex/fs);

    x = zeros(1,L);
    x(indexImpulses) = 1;

    % amplitude modulation
    if strcmp(faultType,'inner')

    if length(fc) > 1
        thetaTime = zeros(1,length(fr));
        for index = 2:length(fr)
            thetaTime(index) = thetaTime(index - 1) + (2*pi/N)/(2*pi*fr(index));
        end
        fcTime = interp1(thetaTime,fc,t,'spline');
        fdTime = interp1(thetaTime,fd,t,'spline');
        fmTime = interp1(thetaTime,fm,t,'spline');

        q = 1 + qAmpMod * cos(2*pi*fcTime.*t + 2*pi*fdTime.*(cumsum(cos(2*pi*fmTime.*t)/fs)));
    else
    q = 1 + qAmpMod * cos(2*pi*fc*t + 2*pi*fd*(cumsum(cos(2*pi*fm*t)/fs)));
    end
    x = q .* x;
    end

    [sdofRespTime] = sdofResponse(fs,k,zita,fn,Lsdof);
    x = fftfilt(sdofRespTime,x);

    L = length(x);
    rng('default'); %set the random generator seed to default (for comparison only)
    SNR = 10^(SNR_dB/10); %SNR to linear scale
    Esym=sum(abs(x).^2)/(L); %Calculate actual symbol energy
    N0 = Esym/SNR; %Find the noise spectral density
    noiseSigma = sqrt(N0); %Standard deviation for AWGN Noise when x is real
    nt = noiseSigma*randn(1,L);%computed noise
xNoise = x + nt; %received signal
 
end