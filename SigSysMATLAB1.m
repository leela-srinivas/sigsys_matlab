%% MATLAB HOMEWORK 1

% Introduction
% * Author: Leela Srinivas
% * Class: ESE 351
% * Date: Created 1/18/2023, Last Edited 1/21/2023

% SET UP INPUT FUNCTIONS

R = 1000;
C = 0.000005;
tau = R * C;

% set up time vector
time_length = 15 * tau;
fs = 44100;
dt = 1/fs;
k = time_length/dt;
k = round(k);

time = linspace(0, time_length, k);

% set up delayed step function input
t0 = 2 * tau; 
u0 = zeros(1, k);

for i = 1:k
    if time(i) < t0
        u0(i) = 0;
    else
        u0(i) = 1;
    end
end

% set up sinusoidal inputs, time vectors for them
A = 5;
w1 = 2 * pi * 10;
time_final_sin1 = 15 * 1/10;
w2 = 2 * pi * 3200;
time_final_sin2 = 15 * 1/3200;
phi1 = 0;
phi2 = pi/2;

j1 = time_final_sin1 * fs;
j2 = cast(time_final_sin2 * fs, "int8");

time_sin1 = linspace(0, time_final_sin1, j1);
time_sin2 = linspace(0, time_final_sin2, j2);

sinusoid1 = A * cos(w1 * time_sin1 - phi1);
sinusoid2 = A * cos(w2 * time_sin2 - phi2); 

% CALCULATING AND PLOTTING OUTPUT VOLTAGES, RESISTIVE LOAD RC CIRCUIT

%set up vectors for output voltages
vout1 = zeros(1, k);
vout2 = zeros(1, j1);
vout3 = zeros(1, j2);

figure()

% plotting output of delayed step function input
subplot(3, 1, 1)
plotting_resistive(u0, vout1, time, k, dt, tau);
title("Calculated Response to Delayed Step Input, with Resistive Load")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")

% plotting output of sinusoid 1 input
subplot(3, 1, 2)
plotting_resistive(sinusoid1, vout2, time_sin1, j1, dt, tau);
title("Calculated Response to Low Frequency AC Source, With Resistive Load")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")

%plotting output of sinusoid 2
subplot(3, 1, 3)
plotting_resistive(sinusoid2, vout3, time_sin2, j2, dt, tau);
title("Calculated Response to High Frequency AC Source, With Resistive Load")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")

% CALCULATING AND PLOTTING OUTPUT VOLTAGES, CAPACITIVE LOAD RC CIRCUIT

vout1_2 = zeros(1, k);
vout2_2 = zeros(1, j1);
vout3_2 = zeros(1, j2);

figure()

% plotting output of delayed step function input
subplot(3, 1, 1)
plotting_capacitive(u0, vout1_2, time, k, dt, tau);
title("Calculated Response to Delayed Step Input Source, With Capacitive Load")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")

% plotting output of sinusoid 1 input
subplot(3, 1, 2)
plotting_capacitive(sinusoid1, vout2_2, time_sin1, j1, dt, tau);
title("Calculated Response to Low Frequency AC Source, With Capacitive Load")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")

% plotting output of sinusoid 2 input
subplot(3, 1, 3)
plotting_capacitive(sinusoid2, vout3_2, time_sin2, j2, dt, tau);
title("Calculated Response to High Frequency AC Source, With Capacitive Load")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")

% VERIFYING WITH LSIM 

% LSIM RESISTIVE LOAD RC CIRCUIT
b = [tau, 0];
a = [tau, 1];
sys = tf(b, a);

figure()

% plotting output of delayed step function input
subplot(3, 1, 1)
lsimming(sys, u0, time);
title("Simulated Resistive Load Response to Delayed Step Input")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")
hold off

% plotting output of sinusoid 1 input
subplot(3, 1, 2)
lsimming(sys, sinusoid1, time_sin1);
title("Simulated Resistive Load Response to Low Frequency Input")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")
hold off

% plotting output of sinusoid 2 input
subplot(3, 1, 3)
lsimming(sys, sinusoid2, time_sin2);
title("Simulated Resistive Load Response to High Frequency Input")
legend("V_{out}", "V_{in}")
ylabel("Voltage (V)")
hold off

% LSIM CAPACITATIVE LOAD RC CIRCUIT
b = 1/tau;
a = [1, 1/tau];
sys = tf(b, a);

figure()

% plotting output of delayed step function input
subplot(3, 1, 1)
lsimming(sys, u0, time);
title("Simulated Capacitive Load Response to Delayed Step Input")
hold off

% plotting output of sinusoid 1 input
subplot(3, 1, 2)
lsimming(sys, sinusoid1, time_sin1);
title("Simulated Capacitive Load Response to Low Frequency Input")
hold off

% plotting output of sinusoid 2 input
subplot(3, 1, 3)
lsimming(sys, sinusoid2, time_sin2);
title("Simulated Capacitive Load Response to High Frequency Input")
hold off

% PROCESSING AUDIO SIGNAL

%loading audio signal from .mat file into array
load supercut.mat
filename = 'supercut.mp3';
Fs = fs;
audiowrite(filename,supercut_sig, Fs);
clear supercut_sig
[supercut_sig,Fs] = audioread('supercut.mp3');

% separating right and left channel of signal (right was used, arbitrarily,
% given they sound identical)
supercut_sig_right = supercut_sig(:, 1);
supercut_sig_left = supercut_sig(:, 2);

numData = length(supercut_sig_right);

% setting up time vector
time_fin_supercut = numData/Fs;
time_supercut = linspace(0, time_fin_supercut, numData);

% plotting right and left channel of signal
figure();
subplot(2, 1, 1);
plot(time_supercut, supercut_sig_right);
title("Right Channel Signal");
xlim([0 time_fin_supercut])
ylabel("Amplitude")
xlabel("Time (s)")
subplot(2, 1, 2);
plot(time_supercut, supercut_sig_left);
title("Left Channel Signal");
ylabel("Amplitude")
xlabel("Time (s)")
xlim([0 time_fin_supercut])

% setting up output signal vectors
right_output_resistive = zeros(1, numData);
right_output_resistive2 = zeros(1, numData);
right_output_capacitive = zeros(1, numData);


% processing signal through resistive load circuit with original time
% constant

R = 1000;
C = 0.000005;
tau = R * C;

figure()
subplot(2, 1, 1)
right_output_resistive = plotting_resistive(supercut_sig_right, right_output_resistive, time_supercut, numData, dt, tau);
title(["Right Channel Signal after Processing via Resistive Load with RC = ", sprintf('%.3f',tau)]);
xlim([0 time_fin_supercut])
legend("Input", "Output")
ylabel("Amplitude")


% processing signal through resistive load circuit with new time
% constant to observe attenuation
R = 5;
C = 0.000005;
tau = R * C;

subplot(2, 1, 2)
right_output_resistive2 = plotting_resistive(supercut_sig_right, right_output_resistive2, time_supercut, numData, dt, tau);
title(["Right Channel Signal after Processing via Resistive Load with RC = ", sprintf('%.6f',tau)]);
xlim([0 time_fin_supercut])
legend("Input", "Output")
ylabel("Amplitude")

% processing signal through capacitative load circuit 
R = 1000;
C = 0.000005;
tau = R * C;

figure()
right_output_capacitive = plotting_capacitive(supercut_sig_right, right_output_capacitive, time_supercut, numData, dt, tau);
xlim([0 time_fin_supercut])
legend("Input", "Output")
ylabel("Amplitude")
title(["Right Channel Signal after Processing via Capacitive Load with RC = ", sprintf('%.3f',tau)]);

% delay of 1.5 seconds
silence = zeros(1, 1.5 * Fs);

% playing all four signals consecutively, with 1.5 second delay in between
four_signals = [supercut_sig_right.', silence, right_output_resistive, silence, right_output_resistive2, silence, right_output_capacitive];
sound(four_signals, Fs)

function output = plotting_resistive(input, output, time, end_time, dt, tau)
    for i = 2:end_time
        output(i) = input(i) - input(i - 1) + (1 - dt/tau)*output(i - 1);
    end
    
    hold on
    
    plot(time, input, '--')
    plot(time, output)
    xlabel("Time (s)")
end

function output = plotting_capacitive(input, output, time, end_time, dt, tau)
    for i = 2:end_time
        output(i) = (1 - dt/tau)*output(i - 1) + (dt/tau) * input(i - 1);
    end
    
    hold on
    
    plot(time, input, '--')
    plot(time, output)
    xlabel("Time (s)")  
end

function y = lsimming(sys, input, time)
    y1 = lsim(sys, input, time);

    hold on
    plot(time, y1)
    plot(time, input, '--')
    legend("V_{out}", "V_{in}")
    ylabel("Voltage (V)")
end
