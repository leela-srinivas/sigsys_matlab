%% MATLAB HOMEWORK 3

% Introduction
% * Author: Leela Srinivas
% * Class: ESE 351
% * Date: Created 2/18/2024, Last Edited 2/21/2024

R = 1000;
C = 0.000005;
fs = 44100;
dt = 1/fs;

% PART 1 
% 10 HZ
f = 10;
endtime = 3 * 1/f; % defining time vector
time = linspace(0, endtime, fs * endtime + 1);

% defining complex exponential input
omega = 2 * pi * f;
exp_input1 = exp(1j * omega * time);

% high pass filter: lsim
tau = R * C;
b = [tau, 0];
a = [tau, 1];
sys_hp = tf(b, a);

% plotting response from high pass filter
y1 = lsim(sys_hp, exp_input1, time);
figure();
subplot(3, 1, 1)
plotting(time, real(exp_input1), real(y1));
title("High Pass RC Circuit Response, f = 10 Hz, Real Component");
subplot(3, 1, 2)
plotting(time, imag(exp_input1), imag(y1));
title("High Pass RC Circuit Response, f = 10 Hz, Imaginary Component");
subplot(3, 1, 3)
plotting(time, abs(exp_input1), abs(y1));
title("High Pass RC Circuit Response, f = 10 Hz, Magnitude");
% Here, the signal has been attenuated, which makes sense because it's a
% high pass filter and the frequency is low. The magnitude in the steady
% state is nearly zero, which makes sense because the signal has been
% attenuated, so there isn't a large output voltage.

% low pass filter: lsim
b = 1/tau;
a = [1, 1/tau];
sys_lp = tf(b, a);

% plotting response from high pass filter
y2 = lsim(sys_lp, exp_input1, time);
figure();
subplot(3, 1, 1)
plotting(time, real(exp_input1), real(y2));
title("Low Pass RC Circuit Response, f = 10 Hz, Real Component");
subplot(3, 1, 2)
plotting(time, imag(exp_input1), imag(y2));
title("Low Pass RC Circuit Response, f = 10 Hz, Imaginary Component");
subplot(3, 1, 3)
plotting(time, abs(exp_input1), abs(y2));
title("Low Pass RC Circuit Response, f = 10 Hz, Magnitude");
% Here, the signal is nearly identical to the input, which makes sense
% because it's a low pass filter and the frequency is low. The steady-state
% magnitude is nearly 1, which makes sense because signal was barely
% attenuated, so the output voltage should be near the input voltage.

% 1000 Hz
f = 1000;
endtime = 20 * 1/f;
time = linspace(0, endtime, fs * endtime + 1);

% defining complex exponential input
omega = 2 * pi * f;
exp_input2 = exp(1j * omega * time);

% high pass
y3 = lsim(sys_hp, exp_input2, time);
figure();
subplot(3, 1, 1)
plotting(time, real(exp_input2), real(y3));
title("High Pass RC Circuit Response, f = 1000 Hz, Real Component");
subplot(3, 1, 2)
plotting(time, imag(exp_input2), imag(y3));
title("High Pass RC Circuit Response, f = 1000 Hz, Imaginary Component");
subplot(3, 1, 3)
plotting(time, abs(exp_input2), abs(y3));
title("High Pass RC Circuit Response, f = 1000 Hz, Magnitude");
% Here, the signal is nearly identical to the input, which makes sense
% because it's a high pass filter and the frequency is high. The steady-state
% magnitude settles to 1, which makes sense because signal was barely
% attenuated, so the output voltage should be near the input voltage.

% low pass filter
y4 = lsim(sys_lp, exp_input2, time);
figure();
subplot(3, 1, 1)
plotting(time, real(exp_input2), real(y4));
title("Low Pass RC Circuit Response, f = 1000 Hz, Real Component");
subplot(3, 1, 2)
plotting(time, imag(exp_input2), imag(y4));
title("Low Pass RC Circuit Response, f = 1000 Hz, Imaginary Component");
subplot(3, 1, 3)
plotting(time, abs(exp_input2), abs(y4));
title("Low Pass RC Circuit Response, f = 1000 Hz, Magnitude");
% Here, the signal has been attenuated, which makes sense because it's a
% low pass filter and the frequency is high. The magnitude in the steady
% state is nearly zero, which makes sense because the signal has been
% attenuated, so there isn't a large output voltage.

% PART 2

% getting frequencies, finding frequency responses of each filter
frequencies = logspace(1, 4, 9); 
H_hp = zeros(1, length(frequencies));
H_lp = zeros(1, length(frequencies));

for i=1:length(frequencies)
    H_hp(i) = part_two_frequency(frequencies(i), sys_hp);
    H_lp(i) = part_two_frequency(frequencies(i), sys_lp);
end

% converting to db
db_hp = mag2db(abs(H_hp));
db_lp = mag2db(abs(H_lp));

% finding the phase, normalizing
phase_normal_hp = angle(H_hp)/pi;
phase_normal_lp = angle(H_lp)/pi;

figure();
subplot(2, 1, 1);
plot(frequencies, db_hp); xscale log
title("Magnitude of Frequency Response, High Pass Filter")
ylabel("Magnitude (db)"); xlabel("Frequency (Hz)");
subplot(2, 1, 2);
plot(frequencies, phase_normal_hp); xscale log
title("Phase of Frequency Response, High Pass Filter")
ylabel("Phase / \pi"); xlabel("Frequency (Hz)");
% the magnitude of the frequency response increases as the frequency
% increases, which makes sense because the high pass filter should preserve
% the input signal's magnitude more when the signal's frequency is higher.
% The phase decreases as the frequency increases, asymptotically approaching 
% zero, which makes sense because the phase is proportional to the arctan 
% of the reciprocal of the frequency, which decreases like so.

figure();
subplot(2, 1, 1);
plot(frequencies, db_lp); xscale log
title("Magnitude of Frequency Response, Low Pass Filter")
ylabel("Magnitude (db)"); xlabel("Frequency (Hz)");
subplot(2, 1, 2);
plot(frequencies, phase_normal_lp); xscale log
title("Phase of Frequency Response, Low Pass Filter")
ylabel("Phase / \pi"); xlabel("Frequency (Hz)");
% the magnitude of the frequency response decreases as the frequency
% increases, which makes sense because the low pass filter attenuates the 
% higher frequency signals, therefore decreasing the magnitude of H(jw).
% The phase increases in magnitude as the frequency increases, asymptotically 
% approaching a negative value, which makes sense because the phase is 
% -arctan(2 * pi * f * R * C), so it'll asymptotically approach a negative
% as the frequency increases.

% PART THREE
fs = 10000;
R = 4800; C = 0.00000025;
f = 60;

endtime = 3;
time = linspace(0, endtime, fs * endtime + 1);

% generating square wave
square_wave = square(2 * pi * f * time);

% low pass filter
tau = R * C;
b = 1/tau;
a = [1, 1/tau];
sys_lp_p3 = tf(b, a);

% high pass filter
b = [tau, 0];
a = [tau, 1];
sys_hp_p3 = tf(b, a);

% found each filter's response to the square wave to find out which makes
% the square wave more like a sinusoid.
square_through_lp = lsim(sys_lp_p3, square_wave, time);
square_through_hp = lsim(sys_hp_p3, square_wave, time);

figure();
hold on
plot(time, square_through_lp);
plot(time, square_through_hp);
plot(time, square_wave);
hold off
xlim([0 6/60]);
legend("Low Pass", "High Pass", "Input");
ylabel("Voltage (V)"); xlabel("Time (s)");
title("RC Circuits' Responses to Square Wave Input");

% low pass filter makes the square wave closer to being sinusoidal, so I
% used a cascade of low pass filters. Here, I experiment with how many low
% pass filters I use in my cascade.
lp_2x = filter_n_times(2, square_wave, sys_lp_p3, time);
lp_3x = filter_n_times(3, square_wave, sys_lp_p3, time);
lp_4x = filter_n_times(4, square_wave, sys_lp_p3, time);

figure();
subplot(3, 1, 1);
hold on
plot(time, square_wave);
plot(time, lp_2x);
hold off
xlim([0 6/60]);
legend("Input", "Output");
ylabel("Voltage (V)"); xlabel("Time (s)");
title("Response to Filtering Square Wave through Low Pass Twice");
subplot(3, 1, 2);
hold on
plot(time, square_wave);
plot(time, lp_3x);
hold off
xlim([0 6/60]);
legend("Input", "Output");
ylabel("Voltage (V)"); xlabel("Time (s)");
title("Response to Filtering Square Wave through Low Pass Thrice");
subplot(3, 1, 3);
hold on
plot(time, square_wave);
plot(time, lp_4x);
hold off
xlim([0 6/60]);
legend("Input", "Output");
ylabel("Voltage (V)"); xlabel("Time (s)");
title("Response to Filtering Square Wave through Low Pass Four Times");

% There's tradeoff between how sinusoidal it is and amplitude. Though
% running it through more filters may make the signal more sinusoidal,
% it'll reduce the magnitude of the frequency response, decreasing the
% power efficiency. I chose n = 3 to compromise. 

% Bode plot of filter's frequency response
frequencies = logspace(1, 4, 9); 
H_p3 = zeros(1, length(frequencies));

for i=1:length(frequencies)
    H_p3(i) = part_three_frequency(frequencies(i), sys_lp_p3);
end

db_p3 = mag2db(abs(H_p3));
phase_normal_p3 = angle(H_p3)/pi;

figure();
subplot(2, 1, 1);
plot(frequencies, db_p3); xscale log
title("Magnitude of Frequency Response, DC/AC Converter")
ylabel("Magnitude (db)"); xlabel("Frequency (Hz)");
subplot(2, 1, 2);
plot(frequencies, phase_normal_p3); xscale log
title("Phase of Frequency Response, DC/AC Converter")
ylabel("Phase / \pi"); xlabel("Frequency (Hz)");
% The frequency increases exponentially in magnitude, then starts to ricochet 
% around -40 dB. It's unclear if the magnitude approaches a certain amount
% as the frequency increases, as higher frequencies return NaN. The phase
% also decreases from zero exponentially, then shoots up and decreases,
% then decreases further. Both seem like they're settling around a certain
% quantity. It's possible that with the cascade of low pass filters, we
% see each low pass filter act at different "stages" of the system. 

P_in = square_wave.^2;
P_out = lp_3x.^2;

length_p3 = length(P_in);
P_avg_in = (1/length_p3) * dot(P_in, ones(1, length_p3));
P_avg_out = (1/length_p3) * dot(P_out, ones(1, length_p3));

efficiency = P_avg_out/P_avg_in;

a_neg_1 = (1/length_p3) * dot((square_wave .* exp(1j * 120 * pi * time)), ones(1, length_p3));
a_pos_1 = (1/length_p3) * dot((square_wave .* exp(-1j * 120 * pi * time)), ones(1, length_p3));

max_power = abs(a_neg_1)^2 + abs(a_pos_1)^2;
theoretical_max_efficiency = max_power/P_avg_in;

text = ['System efficiency: ', num2str(efficiency), '   Theoretical max efficiency: ', num2str(theoretical_max_efficiency)];
disp(text);

r = theoretical_max_efficiency/efficiency;
disp(['Ideal system is ', num2str(r), ' times as efficent as my system.']);
% it makes sense that my system wasn't as ideal as the ideal system because
% I did sacrifice some efficiency for having an output signal that is closer 
% to being an AC signal.

function y = plotting(time, input, output)
    hold on
    plot(time, output)
    plot(time, input, '--')
    hold off
    legend("V_{out}", "V_{in}"); ylabel("Voltage (V)"); xlabel("Time (s)")
    
end

function H = part_two_frequency(freq, sys)
    f = freq;
    fs = 44100;

    endtime = 5 * 1/f;
    time = linspace(0, endtime, fs * endtime + 1);

    % defining complex exponential input
    omega = 2 * pi * f;
    input = exp(1j * omega * time);

    % output
    y = lsim(sys, input, time);
    
    k = length(input);
    H = (y(k)/input(k));

end

function H = part_three_frequency(freq, sys)
    f = freq;
    fs = 10000;

    endtime = 7 * 1/f;
    time = linspace(0, endtime, fs * endtime + 1);

    % defining complex exponential input
    omega = 2 * pi * f;
    input = exp(1j * omega * time);

    k = length(input);
    input_last = input(k);

    % output

    for i = 1:3
        input = lsim(sys, input, time);
    end
    
    H = (input(k)/input_last);

end

function out = filter_n_times(n, input, sys_lp, time)
    
    for i = 1:n
        input = lsim(sys_lp, input, time);
    end

    out = input;
end