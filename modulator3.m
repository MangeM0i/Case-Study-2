%%
% Define parameters
Tp = 0.1;  % Half pulse width
Ts = 2 * Tp;  % Symbol period is twice the pulse width
dt = Tp / 50;  % Sampling period
fs = 1 / dt;  % Sampling frequency
%N = 50;  % Number of bits per signal
%bits = randi([0 1], 3, N) * 2 - 1;  % Generate random bits (-1, 1) for 3 signals


message = 'I love ESE 351! :)';
binary = str2num(reshape(dec2bin(message)',1,[])');
messageOut = char(bin2dec(num2str(reshape(binary,7,[])')))';
% Define carrier frequencies for the three bands
frequencies = 2*[20, 30, 40];  % Carrier frequencies in Hz
binary = transpose([binary,binary,binary]);
bits = binary;
N = length(binary)
% Define the number of symbol periods to span the window
num_periods = 1;  % Increasing this number will widen the pulse Value of 1 means a gaussian approximation for sinc.
sigma = 0.1;  % Noise standard deviation
% Updated time vector
t = linspace(-num_periods * Ts / 2, num_periods * Ts / 2, num_periods * Ts * fs);

% Recompute the windowed sinc pulse
p = sinc(t / Ts) .* hamming(length(t))';

% Compute Fourier Transform of the pulse for bandwidth calculation
P_f = fft(p, 1024);
P_f = fftshift(P_f);  % Shift zero frequency to center
frequency_vector = linspace(-fs/2, fs/2, 1024);

% Find the -3 dB bandwidth
mag_P_f = abs(P_f);
half_max = max(mag_P_f) / sqrt(2);
index_bw = find(mag_P_f >= half_max);
bw_3dB = frequency_vector(index_bw(end)) - frequency_vector(index_bw(1))*2;  % -3 dB bandwidth

fprintf('The -3 dB bandwidth of the pulse is approximately %.2f Hz\n', bw_3dB);

% Generate the signal vector for the entire message for each band
signal_vectors = zeros(3, ceil(N * Ts * fs));  % Initialize signal vectors

% Assign values of bits at the appropriate locations for each signal
for j = 1:3
    for i = 1:N
        start_idx = round((i - 1) * Ts * fs) + 1;
        end_idx = round(i * Ts * fs);
        signal_vectors(j, start_idx:end_idx) = bits(j, i);
    end
end

% Convolve with the pulse shape to get the transmitted signal y(t) for each band
y_t = zeros(3, length(signal_vectors));
for j = 1:3
    y_t(j, :) = conv(signal_vectors(j, :), p, 'same');
    y_t(j, :) = y_t(j, :) / max(abs(y_t(j, :)));  % Normalizing y_t
end

% Up-conversion
t_up = (0:size(y_t, 2) - 1) / fs;  % Time vector for up-converted signal
upconverted_signals = zeros(size(y_t));
for j = 1:3
    upconverted_signals(j, :) = y_t(j, :) .* cos(2 * pi * frequencies(j) * t_up);
end

% Combine signals for transmission
combined_signal = sum(upconverted_signals, 1);

noise = sigma * randn(1, length(combined_signal));  % Generate noise
received_signal = combined_signal + noise;
noisePower = sum(noise.^2) / length(noise);
signalPower = sum(combined_signal.^2) / length(combined_signal);

% Calculate SNR
snr =(signalPower / noisePower);

% Down-conversion and Matched Filtering
decoded_bits = zeros(size(bits));
wrongBitsMF = zeros(1, 3);
for j = 1:3
    % Down-conversion
    downconverted_signal = received_signal .* cos(2 * pi * frequencies(j) * t_up);

    % Matched Filter receiver
    matched_filter_output = conv(downconverted_signal, fliplr(p), 'same');

    % Sampling and decision making
    for i = 1:N
        sample_time = round(i * Ts * fs);
        sampled_value = matched_filter_output(sample_time);
        if sampled_value > 0
            decoded_bits(j, i) = 1;
        else
            decoded_bits(j, i) = 0;
        end
        if decoded_bits(j, i) ~= bits(j, i)
            wrongBitsMF(j) = wrongBitsMF(j) + 1;
        end
    end
end

% Calculate and display error rates
errorRatesMF = wrongBitsMF / N;
fprintf('Matched Filter Error Rates: %.2f%%, %.2f%%, %.2f%%\n', errorRatesMF * 100);
fprintf('SNR: %.2f \n', snr);
messageIn = char(bin2dec(num2str(reshape(decoded_bits,7,[])')))'

%% Plotting:

% Time vector for plotting
t_plot = (0:length(combined_signal)-1) / fs;

% Plot pulse shapes
figure;
plot(t, p);
title('Windowed Sinc Pulse Shape');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot individual up-converted signals
figure;
for j = 1:3
    subplot(3,1,j);
    plot(t_plot, upconverted_signals(j, :));
    title(['Up-Converted Signal at ' num2str(frequencies(j)) ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

% Plot combined transmitted signal
figure;
plot(t_plot, combined_signal);
title('Combined Transmitted Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot received signal with noise
figure;
plot(t_plot, received_signal);  
title('Received Signal with noise');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot down-converted and filtered signals
figure;
for j = 1:3
    % Down-convert signal
    downconverted_signal = received_signal .* cos(2 * pi * frequencies(j) * t_plot);
    % Matched Filter receiver
    matched_filter_output = conv(downconverted_signal, fliplr(p), 'same');

    subplot(3,1,j);
    plot(t_plot(1:length(matched_filter_output)), matched_filter_output);
    title(['Matched Filter Output for ' num2str(frequencies(j)) ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

% Plot the error rates for each channel
figure;
bar(errorRatesMF * 100);
title('Error Rates for Each Frequency Band');
xlabel('Frequency Band');
ylabel('Error Rate (%)');
xticklabels({'20 Hz', '30 Hz', '40 Hz'});
grid on;


% Define time vector for bit sampling
sample_times = round((1:N) * Ts * fs);

% Plot sent and received bits for each channel
figure;
for j = 1:3
    subplot(3, 1, j);
    stem(sample_times/fs, bits(j, :), 'filled', 'MarkerSize', 4, 'LineStyle', 'none');
    hold on;
    stem(sample_times/fs, decoded_bits(j, :), 'r*', 'MarkerSize', 4);
    hold off;
    title(['Channel ' num2str(frequencies(j)) ' Hz: Sent (blue) vs. Received (red) Bits']);
    xlabel('Time (s)');
    ylabel('Bit Value');
    ylim([-1.5 1.5]);
    grid on;
    legend('Sent Bits', 'Decoded Bits', 'Location', 'best');
end

sgtitle('Comparison of Sent and Received Bits for Each Channel');
