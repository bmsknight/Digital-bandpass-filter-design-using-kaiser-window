close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%

% Digital Filter Design
%B.Madhushan

%parameters
Ap  = 0.13;         % Max passband ripple
Aa  = 56;           % Min stopband ripple
wp1 = 300;          % Lower passband edge
wp2 = 600;          % Upper passband edge
wa1 = 200;          % Lower stopband edge
wa2 = 750;          % Upper stopband edge
ws  = 2000;         % Sampling frequency

%Bt,wc1,wc2
Bt = min((wp1-wa1),(wa2-wp2));
wc1 = wp1 - Bt/2;
wc2 = wp2 + Bt/2;
T = 2*pi/ws;

%delta
deltap = (10^(Ap/20)-1)/(10^(Ap/20)+1);
deltaa = 10^(-Aa/20);

delta = min(deltap,deltaa);

actualAa = -20*log10(delta);

% alpha
if (actualAa<=21)
    alpha = 0;
elseif (actualAa<=50)
    alpha = 0.5842*(actualAa-21)^0.4 + 0.07886*(actualAa-21);
else
    alpha = 0.1102*(actualAa-8.7);
end

% parameter D
if (actualAa<=21)
    D = 0.9222;
else
    D = (actualAa-7.95)/14.36;
end

%N
Bt = min((wp1-wa1),(wa2-wp2));

N = ceil(ws*D/Bt+1);
if(mod(N,2)==0)
    N=N+1;
end

%IoAlpha
 IoAlpha = 0;
 IoValue =1;
 k =1;
 while (IoValue>1*10^(-6))
    IoAlpha =IoAlpha + IoValue;
    IoValue =((1/factorial(k))*((alpha/2)^k))^2;
    k = k+1;
 end
 
 kaiserTemp = zeros(1,(N+1)/2);
 
 for n = 0:(N-1)/2
    beta = alpha*sqrt(1-((2*n)/(N-1))^2);
     
    %IoBeta
    IoBeta = 0;
    IoValue =1;
    k =1;
    
    while (IoValue>1*10^(-6))
       IoBeta =IoBeta + IoValue;
       IoValue =((1/factorial(k))*((beta/2)^k))^2;
       k = k+1;
    end
    kaiserTemp(1,n+1)= IoBeta/IoAlpha;
 end
 
 kaiserWindow = [flip(kaiserTemp),kaiserTemp(2:length(kaiserTemp))];
  
%hnT
hnT = zeros(1,N);
for n = -34:34
    if n==0
        hnT(n+35)=(2/ws)*(wc2-wc1);
    else
        hnT(n+35)=(1/(n*pi))*(sin(wc2*n*T)-sin(wc1*n*T));
    end
end

dFilter = hnT.*kaiserWindow;

%impulse Response
figure;
stem(dFilter,'filled');
title('Impulse Response of the Digital Filter');
xlabel('Samples');
ylabel('Amplitude');

% Magnitude Response
[h,w]=freqz(dFilter);
magnitude = 20*log10(abs(h));
analogFreq = w*ws/(2*pi);

figure;
plot(analogFreq,magnitude);
grid ON;
title('Magnitude Response of the Digital Filter');
xlabel('Ananlog Frequency(rad/s)');
ylabel('Magnitude(dB)');

%Magnitude Response of Passband Ripple
figure;
plot(analogFreq,magnitude);
axis([wp1 wp2 -0.05 0.05]);
grid ON;
title('Magnitude Response of the Pass Band');
xlabel('Ananlog Frequency(rad/s)');
ylabel('Magnitude(dB)');

%Magnitude Response of Lower StopBand
figure;
plot(analogFreq,magnitude);
axis([0 wa1 -120 -50]);
grid ON;
title('Magnitude Response of the Lower Stop Band');
xlabel('Ananlog Frequency(rad/s)');
ylabel('Magnitude(dB)');

%Magnitude Response of Upper StopBand
figure;
plot(analogFreq,magnitude);
axis([wa2 ws/2 -120 -50]);
grid ON;
title('Magnitude Response of the Upper Stop Band');
xlabel('Ananlog Frequency(rad/s)');
ylabel('Magnitude(dB)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Testing

%Input Signal Generation
noOfSamples = 300;
w1 = wa1/2;
w2 = (wp1+wp2)/2;
w3 = (wa2+ws/2)/2;

xnT = zeros(1,noOfSamples);

for i = 1:noOfSamples
    xnT(i) = sin(w1*i*T) + sin(w2*i*T) + sin(w3*i*T);
end

figure;
subplot(3,1,1);
stem(xnT);
title('Input Signal');
xlabel('Samples');
ylabel('amplitude');

%Expected Output Signal
enT = zeros(1,noOfSamples);

for i = 1:noOfSamples
    enT(i) = sin(w2*i*T);
end

subplot(3,1,2);
stem(enT);
axis([0 noOfSamples -2 2]);
title('Expected Output Signal');
xlabel('Samples');
ylabel('amplitude');

%Output Signal

lenfft = length(xnT) + length(dFilter) - 1;

XF = fft(xnT,lenfft);
HF = fft(dFilter,lenfft);
YF = HF.*XF;

ynT = ifft(YF,lenfft);

subplot(3,1,3);
stem(ynT);
axis([0 noOfSamples -2 2]);
title('Output Signal');
xlabel('Samples');
ylabel('amplitude');

%Frequency Response Analysis
w = ws*(-lenfft/2+1:lenfft/2)/lenfft;

XF1 = abs(fftshift(XF));
figure;
subplot(3,1,1);
plot(w,XF1);
axis([-ws/2 ws/2 0 200]);
title('Frequency Spectrum of the Input Signal');
xlabel('Ananlog Frequency(rad/s)');
ylabel('Amplitude');

EF = fft(enT,lenfft);
EF1 = abs(fftshift(EF));
subplot(3,1,2);
plot(w,EF1);
axis([-ws/2 ws/2 0 200]);
title('Frequency Spectrum of the Expected Output Signal');
xlabel('Ananlog Frequency(rad/s)');
ylabel('Amplitude');

YF1 = abs(fftshift(YF));
subplot(3,1,3);
plot(w,YF1);
axis([-ws/2 ws/2 0 200]);
title('Frequency Spectrum of the Output Signal');
xlabel('Ananlog Frequency(rad/s)');
ylabel('Amplitude');
