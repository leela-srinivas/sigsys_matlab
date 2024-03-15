%% MATLAB HOMEWORK 2

% Introduction
% * Author: Leela Srinivas
% * Class: ESE 351
% * Date: Created 1/28/2023, Last Edited 2/2/2023

% initialize constants, any circuit parameters
tau = 0.005;
R = 1000;
C = 0.000005;
fs = 44100;
dt = 1/fs;
endtime = 15 * tau;

% TASK 1a. convolve rectangular function w/ themselves, plot result
% get, convolve, plot rectangle functions with themselves
rect_3 = get_rect(3); 
rect_10 = get_rect(10); 
rect_21 = get_rect(21); 

% convolve rectangle functions w/ themselves
conv1 = conv(rect_3, rect_3); 
conv2 = conv(rect_10, rect_10); 
conv3 = conv(rect_21, rect_21); 

x1 = linspace(0, 43, 44);

% plotting convolved rect. functions
figure();
hold on
stem(x1(1:7), conv1); stem(x1(1:21), conv2); stem(x1(1:43), conv3);
xlabel("n");
ylabel("x[n]")
title("Convolution of Rectangular Function with Themselves");
hold off
legend("N = 3", "N = 10", "N = 21");

% TASK 1b. Get impulse trains, convolve them with last convolution

xaxis = linspace(-200, 200, 401); % domain

% getting impulse train
d1 = get_delta(10); 
d2 = get_delta(25); 
d3 = get_delta(50); 

% plotting impulse trains
figure();
subplot(3, 1, 1); stem(xaxis, d1); 
title("Impulse Train where M = 10"); xlabel("n"); ylabel("Amplitude")
subplot(3, 1, 2); stem(xaxis, d2); 
title("Impulse Train where M = 25") ; xlabel("n"); ylabel("Amplitude")
subplot(3, 1, 3); stem(xaxis, d3); 
title("Impulse Train where M = 50") ; xlabel("n"); ylabel("Amplitude")

% convolve impulse train with x[n] where N = 21
conv4 = conv(conv3, d1(1:362)); 
conv5 = conv(conv3, d2(1:362)); 
conv6 = conv(conv3, d3(1:362)); 

% plotting convolutions with x[n]
figure();
subplot(3, 1, 1); stem(xaxis, conv4(1:401)); 
title("Impulse Train Convolved with Convolution of Rect. Function (M = 10, N = 21)");
xlabel("n"); ylabel("Amplitude");
subplot(3, 1, 2); stem(xaxis, conv5(1:401)); 
title("M = 25"); xlabel("n"); ylabel("Amplitude");
title("Impulse Train Convolved with Convolution of Rect. Function (M = 25, N = 21)");
subplot(3, 1, 3); stem(xaxis, conv6(1:401)); 
title("M = 50"); xlabel("n"); ylabel("Amplitude");
title("Impulse Train Convolved with Convolution of Rect. Function (M = 50, N = 21)");

% TASK 2a. Use filter to simulate DT systems


% set up any constants, time vectors
length_time = round(endtime * fs);
u0 = zeros(1, length_time);
t0 = 2 * tau * fs;

time = linspace(0, endtime, fs * endtime + 1);

% making step input
for i = 1:length_time
    if i < t0
        u0(i) = 0;
    else
        u0(i) = 1;
    end
end

% Making delta function
delta_b = zeros(1, length_time);
delta_b(1) = 1;

% Use filter() to find impulse response for each circuit
% high pass
b_hp_filter = [1, -1]; a_hp_filter = [1, -(1-dt/tau)];
y_hp_dt_impulse = filter(b_hp_filter, a_hp_filter, delta_b); 
% low pass
b_lp_filter = dt/(tau); a_lp_filter = [1 -(1-dt/tau)];
y_lp_dt_impulse = filter(b_lp_filter, a_lp_filter, delta_b); 

% impulse response in continuous time
h_lowpass = (1/tau) * exp(-time/tau) * 1;
h_highpass = - (1/tau) * exp(-time/tau) * 1;
h_highpass(1) = 1 - (1/tau) * exp(-time(1)/tau) * 1;

% plotting found impulse response with given impulse response
% The high pass filter's response isn't as visible, it's zoomed in so we
% can see it better, given response is scaled down by 200
figure();
subplot(2, 2, 1);
stem(y_hp_dt_impulse); 
title("Impulse Response of High Pass RC Filter, Calculated with filter()")
xlabel("n"); ylabel("Amplitude")
subplot(2, 2, 2);
stem(y_hp_dt_impulse);
xlim([0 300]); ylim([-0.01, 0.01]);
title("Impulse Response of High Pass RC Filter, Calculated with filter() and Zoomed In")
xlabel("n"); ylabel("Amplitude")
subplot(2, 2, [3 4]);
plot(time, (h_highpass) * 1/200);
xlabel("Time (s)"); ylabel("Amplitude")
title("Impulse Response of High Pass RC Filter, Given Equation")

% same for low pass filter, given response scaled down by 200
figure();
subplot(2, 1, 1);
stem(y_lp_dt_impulse); 
title("Impulse Response of Low Pass RC Filter, Calculated with filter()")
xlabel("n"); ylabel("Amplitude")
subplot(2, 1, 2);
plot(time, (h_lowpass) * 1/200);
xlabel("Time (s)"); ylabel("Amplitude")
title("Impulse Response of Low Pass RC Filter, Given Equation")

% TASK 2B: compute response of each circuit to shifted step function input
% i. compute via filter
y_hp_dt_step = filter(b_hp_filter, a_hp_filter, u0);
y_lp_dt_step = filter(b_lp_filter, a_lp_filter, u0);

% ii. compute via convolution with impulse response found in step 2a.
y_lowpass = conv(u0, y_lp_dt_impulse); y_highpass = conv(u0, y_hp_dt_impulse);

% iii. compute via lsim
% high pass
b_hp_lsim = [tau, 0];
a_hp_lsim = [tau, 1];
sys1 = tf(b_hp_lsim, a_hp_lsim);
y_highpass_lsim = lsim(sys1, u0, time);

% low pass
b = 1/tau;
a = [1, 1/tau];
sys2 = tf(b, a);
y_lowpass_lsim = lsim(sys2, u0, time);

% plot lsim result, conv result, filter() result
figure();

% high pass
subplot(2, 2, 1)
hold on
plot(y_hp_dt_step); plot(y_highpass(1:length_time)); 
hold off
title("Highpass Filter Step Response, DT"); 
legend("Calculated with filter()", "Calculated with Convolution");
xlabel("n"); ylabel("y[n]");
subplot(2, 2, 2); plot(time, y_highpass_lsim);
title("High Pass Filter Step Response, Found with lsim, CT"); 
xlabel("Time (s)"); ylabel("y(t)");
% low pass
subplot(2, 2, 3)
hold on
plot(y_lp_dt_step); plot(y_lowpass(1:length_time));
title("Lowpass Filter Step Response, DT"); 
legend("Calculated with filter()", "Calculated with Convolution")
xlabel("n"); ylabel("y[n]");
hold off
subplot(2, 2, 4);
plot(time, y_lowpass_lsim);
title("Lowpass Filter Step Response, Found with lsim, CT");
xlabel("Time (s)"); ylabel("y(t)");

% TASK 2c. Find rect. functions' response to step  
% concatenate zero vector with rect function for filter()
rect_3_zeros = [rect_3, zeros(1, 4)];
rect_10_zeros = [rect_10, zeros(1, 11)];
rect_21_zeros = [rect_21, zeros(1, 22)];

% like convolving rectangular pulse with itself
y_conv3 = filter(rect_3_zeros, 1, rect_3_zeros); 
y_conv10 = filter(rect_10_zeros, 1, rect_10_zeros);
y_conv21 = filter(rect_21_zeros, 1, rect_21_zeros); 

% convolve result with impulse train to get response
y_21_d1 = conv(y_conv21, d1(1:361)); 
y_21_d2 = conv(y_conv21, d2(1:361));
y_21_d3 = conv(y_conv21, d3(1:361));

% plotting results
figure(); 
subplot(3, 1, 1); 
stem(xaxis, y_21_d1(1:401)); 
xlabel("n"); ylabel("Amplitude");
title("Convolution of Rect. Function with Itself, Impulse Train (M = 10)")
subplot(3, 1, 2); 
stem(xaxis, y_21_d2(1:401)); 
title("Convolution of Rect. Function with Itself, Impulse Train (M = 25)")
xlabel("n"); ylabel("Amplitude");
subplot(3, 1, 3);
stem(xaxis, y_21_d3(1:401)); 
xlabel("n"); ylabel("Amplitude");
title("Convolution of Rect. Function with Itself, Impulse Train (M = 50)")

% TASK 3a. Passing chirp through systems

% This code provided by Ethan Cuka and Dr. Jason Trobaugh via 
% https://wustl.instructure.com/courses/121863/assignments/632432?module_item_id=2117432
t = 0:dt:3; % time vector
fmin = 10; fmax = 280; % the corner frequency was found to be 31.8310 Hz, so I reduced fmax considerably
fchirp = (fmax-fmin).*t/max(t)+fmin; % chirp instantaneous frequency
xchirp = cos(2*pi*fchirp/2.*t); % chirp signal

% passing chirp through filters
ychirpHPF = filter(b_hp_filter, a_hp_filter,xchirp);
ychirpLPF = filter(b_lp_filter, a_lp_filter, xchirp);

% plotting original chirp and filtered chirps
figure();
subplot(3,1,1), plot(fchirp, xchirp); title('Chirp Input'); 
xlabel('Frequency (Hz)'); ylabel("Amplitude"); xlim([fmin 110]);
subplot(3,1,2), plot(fchirp,(ychirpHPF)); title('Chirp Filtered via DT Highpass Filter'),
xlabel('Frequency (Hz)'); ylabel("Amplitude"); xlim([fmin 110]);
subplot(3,1,3), plot(fchirp,(ychirpLPF)); title('Chirp Filtered via DT Lowpass Filter'),
xlabel('Frequency (Hz)'); ylabel("Amplitude"); xlim([fmin 110]);

% TASK 3b. Finding DT rectangular functions' response to the chirp
% get rectangular functions
rect_10 = get_rect(10);
rect_25 = get_rect(25);
rect_50 = get_rect(50);
rect_100 = get_rect(100);

% finding frequency response
ychirp10 = filter(rect_10, 1, xchirp);
ychirp25 = filter(rect_25, 1, xchirp);
ychirp50 = filter(rect_50, 1, xchirp);
ychirp100 = filter(rect_100, 1, xchirp);

% plotting filtered chirps
figure();
subplot(4, 1, 1)
plot(fchirp,ychirp10); 
xlabel('Frequency (Hz)'); ylabel("Amplitude"); 
title("Chirp Signal Through Rectangular Function (N = 10)"); xlim([60 fmax]);
subplot(4, 1, 2)
plot(fchirp,ychirp25); xlabel('Frequency (Hz)'); ylabel("Amplitude");
title("Chirp Signal Through Rectangular Function (N = 25)"); xlim([60 fmax]);
subplot(4, 1, 3)
plot(fchirp,ychirp50); xlabel('Frequency (Hz)'); ylabel("Amplitude");
title("Chirp Signal Through Rectangular Function (N = 50)"); xlim([60 fmax]);
subplot(4, 1, 4)
plot(fchirp,ychirp100); xlabel('Frequency (Hz)'); ylabel("Amplitude");
title("Chirp Signal Through Rectangular Function (N = 100)"); xlim([60 fmax]);

% Part 4. Filter audio signal through RC circuits and DT rectangular
% functions
% read supercut signal in
[supercut_sig,Fs] = audioread('supercut.mp3');
supercut_sig_right = supercut_sig(:, 1);
numData = length(supercut_sig_right);

% finding time vector
time_fin_supercut = numData/Fs;
time_supercut = linspace(0, time_fin_supercut, numData);

% pass through RC circuit filters
supercutHPF = filter(b_hp_filter, a_hp_filter,supercut_sig_right);
supercutLPF = filter(b_lp_filter, a_lp_filter, supercut_sig_right);

% pass through high pass filter with different time constant
tau = 0.0001;
b_hp_filter_2 = [1, -1]; a_hp_filter_2 = [1, -(1-dt/tau)];
supercutHPF_2 = filter(b_hp_filter_2, a_hp_filter_2,supercut_sig_right);

% plotting "Supercut" through high pass filter with different time constants
figure();
subplot(2, 1, 1);
hold on
plot(time_supercut, supercut_sig_right);
plot(time_supercut, supercutHPF);
legend("Original Signal", "Filtered Signal")
title("Chosen Signal with High Pass Filter, RC = 0.005");
xlabel("Time(s)"); ylabel("Amplitude");
hold off
subplot(2, 1, 2);
hold on
plot(time_supercut, supercut_sig_right);
plot(time_supercut, supercutHPF_2);
legend("Original Signal", "Filtered Signal")
xlabel("Time(s)"); ylabel("Amplitude");
title("Chosen Signal with High Pass Filter, RC = 0.0001");
hold off

% plotting "Supercut" through low pass filter
figure();
hold on
plot(time_supercut, supercut_sig_right);
plot(time_supercut, supercutLPF);
legend("Original Signal", "Filtered Signal")
xlabel("Time(s)"); ylabel("Amplitude");
title("Chosen Signal with Low Pass Filter, RC = 0.005");
hold off

% Finding response of Rectangular Functions to "Supercut"
supercut10 = filter(rect_10, 1, supercut_sig_right);
supercut25 = filter(rect_25, 1, supercut_sig_right);
supercut50 = filter(rect_50, 1, supercut_sig_right);
supercut100 = filter(rect_100, 1, supercut_sig_right);

% plotting results
figure();
subplot(4, 1, 1)
hold on
plot(time_supercut, supercut10);
plot(time_supercut, supercut_sig_right)
hold off
title("Rectangular Function (N = 10) Response to Audio Signal");
xlabel("Time (s)"); ylabel("Amplitude");
legend("Frequency Response", "Original Signal");
subplot(4, 1, 2)
hold on
plot(time_supercut, supercut25)
plot(time_supercut, supercut_sig_right);
hold off
title("Rectangular Function (N = 25) Response to Audio Signal");
xlabel("Time (s)"); ylabel("Amplitude");
legend("Frequency Response", "Original Signal");
subplot(4, 1, 3)
hold on
plot(time_supercut, supercut50)
plot(time_supercut, supercut_sig_right);
hold off
title("Rectangular Function (N = 50) Response to Audio Signal");
xlabel("Time (s)"); ylabel("Amplitude");
legend("Frequency Response", "Original Signal");
subplot(4, 1, 4)
hold on
plot(time_supercut, supercut100)
plot(time_supercut, supercut_sig_right);
hold off
title("Rectangular Function (N = 100) Response to Audio Signal");
xlabel("Time (s)"); ylabel("Amplitude");
legend("Frequency Response", "Original Signal");

% 1.5 seconds of silence
silence = zeros(1, 1.5 * fs);

% used to hear sound
%sound([supercutHPF.', silence, supercutHPF_2.', silence, supercutLPF.'], Fs);
%sound([supercut10.', silence, supercut25.', silence, supercut50.', silence, supercut100.'], Fs);

snapnow

function rect = get_rect(N)

    rect = zeros(1, N + 1);

    for i=1:N+1
        if i <= N + 1
            rect(i) = 1;
        else
            rect(i) = 0;
        end
    end
end

function delta = get_delta(M)

    delta = zeros(1, 401);
    
    for k = 1:401
        if mod(k - 1, M) == 0
            delta(k) = 1;
        else
            delta(k) = 0;
        end
    end
end