%% Define half pulse width (Tp) and sampling period (dt)
Tp = 0.1;
dt = Tp/50;

% Define time vector for pulse
t = -Tp:dt:Tp;

% Define triangular pulse function
p_t = (1 - abs(t)/Tp).*(abs(t) < Tp);

% Numerical approximation of Fourier Transform using FFT
P_f = fft(p_t); 
fs= 100;
f = linspace(-fs/2,fs/2,length(P_f)); % frequency vector for FFT output

%Renamed N to bitsT
bitsT=50;
bits = 2*((rand(1,bitsT)>0.5)-0.5);

% Define bit rate (fb) based on our choice (1/Tp or 1/(2*Tp))
fb = 1/Tp; 

% Symbol period (Ts) based on bit rate
Ts = 1/fb;
t2 = 0:dt:(bitsT*Ts);

% Base zero vector for bit insertion
bSpaced = zeros(1, length(t2)); 

% Loop for bit rate adjusted bit vector
i = 1;
j = 1;
for loop_index = 0:length(bits)-1 

    bSpaced(i) = bits(j); 
    i = round(Ts/dt)*j;    
    j = j + 1;   
end

y_t=conv(bSpaced,p_t);

tPlot=0:dt:dt*(length(y_t)-1);

% Noise standard deviation 
sigma = 1;
noise = sigma*randn(1,length(y_t));


% Received signal with noise
r_t = y_t + noise;

% Sign-based receiver
signX=zeros(1,length(t2));
i=1;
j=1;
for loop_index = 0:length(bits)-1

if r_t(i)>0
signX(i)=1;

else 
signX(i)=-1;
    
end
i=j*round(Ts/dt);
j=j+1;
end

%Matched filter

z_t=conv(r_t,p_t);
matchedX=zeros(1,length(t2));
i=1;
j=1;
for loop_index = 0:length(bits)-1

if z_t(i)>0
matchedX(i)=1;

else 
matchedX(i)=-1;
    
end
i=j*round(Ts/dt);
j=j+1;
end


Now We Plot our Data

figure;
subplot(2,1,1)
plot(t,p_t)
xlabel('Time (s)')
ylabel('Pulse Amplitude')
title('Triangular Pulse Shape')

subplot(2,1,2)
semilogy(f,abs(P_f)) % we use semilogy to better view the large data scale on the y-axis
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of Triangular Pulse (FFT)')

figure;
subplot(2,1,1)
plot(tPlot,y_t)
xlabel('Time (s)')
ylabel('Signal Amplitude')
ylim([-2 2])
title('Noise Free Transmitted Signal')

subplot(2,1,2)
plot(tPlot,r_t)
xlabel('Time(s)')
ylabel('Signal Amplitude')
title('Noisy Received Signal')

figure;
subplot(2,1,1)
plot(t2,bSpaced)
xlabel('Time (s)')
ylabel('Signal Amplitude')
title('Sent Message Sign Based Receiver')

subplot(2,1,2)
plot(t2,signX)
xlabel('Time (s)')
ylabel('Signal Amplitude')
title('Decoded Message Sign Based Receiver')

figure;
subplot(2,1,1)
plot(t2,bSpaced)
xlabel('Time(s)')
ylabel('Signal Amplitude')
title('Sent Message Matched Filter')

subplot(2,1,2)
plot(t2,matchedX)
xlabel('Time(s)')
ylabel('Signal Amplitude')
title('Decoded Message Matched Filter')

SNR = sum(y_t.^2) / sum(noise.^2);
ER_sign = (sum(abs(bSpaced - signX)) / length(bSpaced));
ER_Matched = (sum(abs(bSpaced - matchedX)) / length(bSpaced));

fprintf('Bit Rate (Hz): %f\n', fb);
fprintf('Noise Std. Dev. (sigma): %f\n', sigma);
fprintf('SNR (dB): %.2f\n', SNR);
fprintf('Error Rate (Sign-based): %f\n', ER_sign);
fprintf('Error Rate (Matched Filter): %f\n', ER_Matched);


Task 1, bit rate = 1/(2Tp)
%% Define half pulse width (Tp) and sampling period (dt)
Tp = 0.1;
dt = Tp/50;

% Define time vector for pulse
t = -Tp:dt:Tp;

% Define triangular pulse function
p_t = (1 - abs(t)/Tp).*(abs(t) < Tp);

% Numerical approximation of Fourier Transform using FFT
P_f = fft(p_t); 


%Renamed N to bitsT
bitsT=50;
bits = 2*((rand(1,bitsT)>0.5)-0.5);


% Define bit rate (fb) based on our choice (1/Tp or 1/(2*Tp))
fb = 1/(2*Tp); 

% Symbol period (Ts) based on bit rate
Ts = 1/fb;
t2 = 0:dt:(bitsT*Ts);
% Base zero vector for bit insertion
bSpaced = zeros(1, length(t2)); 

% Loop for bit rate adjusted bit vector
i = 1;
j = 1;
for loop_index = 0:length(bits)-1 

    bSpaced(i) = bits(j); 
    i = round(Ts/dt)*j;    
    j = j + 1;   
end

y_t=conv(bSpaced,p_t);

tPlot=0:dt:dt*(length(y_t)-1);
vec=randn(1,length(y_t));
noise = sigma*vec;

% Received signal with noise
r_t = y_t + noise;

% Sign-based receiver
sign2X=zeros(1,length(t2));
i=1;
j=1;
for loop_index = 0:length(bits)-1

if r_t(i)>0
sign2X(i)=1;

else 
sign2X(i)=-1;
    
end
i=j*round(Ts/dt);
j=j+1;
end

%Matched filter

z_t=conv(r_t,p_t);
matched2X=zeros(1,length(t2));
i=1;
j=1;
for loop_index = 0:length(bits)-1

if z_t(i)>0
matched2X(i)=1;

else 
matched2X(i)=-1;
    
end
i=j*round(Ts/dt);
j=j+1;
end


Now We Plot our Data

figure;
subplot(2,1,1)
plot(t,p_t)
xlabel('Time (s)')
ylabel('Pulse Amplitude')
title('Triangular Pulse Shape')

subplot(2,1,2)
semilogy(f,abs(P_f)) 
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of Triangular Pulse (FFT)')

figure;
subplot(2,1,1)
plot(tPlot,y_t)
xlabel('Time (s)')
ylabel('Signal Amplitude')
ylim([-2 2])
title('Noise Free Transmitted Signal')

subplot(2,1,2)
plot(tPlot,r_t)
xlabel('Time(s)')
ylabel('Signal Amplitude')
title('Noisy Received Signal')

figure;
subplot(2,1,1)
plot(t2,bSpaced)
xlabel('Time (s)')
ylabel('Signal Amplitude')
title('Sent Message Sign Based Receiver')

subplot(2,1,2)
plot(t2,sign2X)
xlabel('Time (s)')
ylabel('Signal Amplitude')
title('Decoded Message Sign Based Receiver')

figure;
subplot(2,1,1)
plot(t2,bSpaced)
xlabel('Time(s)')
ylabel('Signal Amplitude')
title('Sent Message Matched Filter')

subplot(2,1,2)
plot(t2,matched2X)
xlabel('Time(s)')
ylabel('Signal Amplitude')
title('Decoded Message Matched Filter')

SNR2 = sum(y_t.^2) / sum(noise.^2);

%% The error rates shown below  are very small, so it's difficult to see any strong differences. These error rates can be amplified to better visualize the difference between them. The matched filter usually has a smaller error rate.

ER_sign2 = (sum(abs(bSpaced - sign2X)) / length(bSpaced));
ER_Matched2 = (sum(abs(bSpaced - matched2X)) / length(bSpaced));

fprintf('Bit Rate (Hz): %f\n', fb);
fprintf('Noise Std. Dev. (sigma): %f\n', sigma);
fprintf('SNR (dB): %.2f\n', SNR2);
fprintf('Error Rate (Sign-based): %f\n', ER_sign2);
fprintf('Error Rate (Matched Filter): %f\n', ER_Matched2);


Task 2: Performance Analysis
%% Define half pulse width (Tp) and sampling period (dt)
Tp = 0.1;
dt = Tp/50;

% Define time vector for pulse
t = -Tp:dt:Tp;

% Define triangular pulse function
p_t = (1 - abs(t)/Tp).*(abs(t) < Tp);

% Numerical approximation of Fourier Transform using FFT
P_f = fft(p_t); 


%Renamed N to bitsT and bits made constant to analyze sigma change. If this
%is to be tested again, a new constant bits1 value will need to be made.
%Some of the code below is commented out because the constant vectors no
%longer exist.
bitsT=50;
%bits = bits1

% Define bit rate (fb) based on our choice (1/Tp or 1/(2*Tp))
fb = 1/Tp; 

% Symbol period (Ts) based on bit rate
Ts = 1/fb;
t2 = 0:dt:(bitsT*Ts);

% Base zero vector for bit insertion
bSpaced = zeros(1, length(t2)); 

% Loop for bit rate adjusted bit vector
i = 1;
j = 1;
for loop_index = 0:length(bits)-1 

    bSpaced(i) = bits(j); 
    i = round(Ts/dt)*j;    
    j = j + 1;   
end

y_t=conv(bSpaced,p_t);

tPlot=0:dt:dt*(length(y_t)-1);

% Noise standard deviation and random vector made constant to analyze sigma
% change. If this is to be tested again, a new constant vec1 value will need to be made
sigma = .0001;
vec=randn(1,length(y_t));
%vec1
%noise = sigma*vec1;


% Received signal with noise
%r_t = y_t + noise;

% Sign-based receiver
signX=zeros(1,length(t2));
i=1;
j=1;
for loop_index = 0:length(bits)-1

if r_t(i)>0
signX(i)=1;

else 
signX(i)=-1;
    
end
i=j*round(Ts/dt);
j=j+1;
end

%Matched filter

z_t=conv(r_t,p_t);
matchedX=zeros(1,length(t2));
i=1;
j=1;
for loop_index = 0:length(bits)-1

if z_t(i)>0
matchedX(i)=1;



else 
matchedX(i)=-1;
    
end
i=j*round(Ts/dt);
j=j+1;
end

SNR = sum(y_t.^2) / sum(noise.^2);
ER_sign = (sum(abs(bSpaced - signX)) / length(bSpaced));
ER_Matched = (sum(abs(bSpaced - matchedX)) / length(bSpaced));

fprintf('Bit Rate (Hz): %f\n', fb);
fprintf('Noise Std. Dev. (sigma): %f\n', sigma);
fprintf('SNR (dB): %.2f\n', SNR);
fprintf('Error Rate (Sign-based): %f\n', ER_sign);
fprintf('Error Rate (Matched Filter): %f\n', ER_Matched);
