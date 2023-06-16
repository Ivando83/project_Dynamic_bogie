%%Program Eksitasi Harmonik
clear all;
clc;

%% Sistem Parameter Bogie
K1    =580000; %Newton
K2    =580000; %Newton
K3    =580000; %Newton
K4    =580000; %Newton
m     =29892;  %Newton
J     =2750;   %Newton
C1    =34000;  %Newton
C2    =34000;  %Newton

%Panjang titik pusat ke titik pusat roda
L1    = 1.5; %Meter
L2    = 1.5; %Meter

%%Amplitude base eksitasi 
amp   = 0.008; %[m]

%% Kecepatan pada Bogie
v     = 16;  % [m/s] dari 60 km/h

%% Panjang eksitasi harmonik antar bantalan
p_exc = 0.55; % [m]

%%eksitasi harmonik amp*cos(omega*time)
omega = 2*pi*v/p_exc; %[rad/s]

%sudut fase antara eksitasi roda depan dan belakang
delta = ((L1 + L2)/p_exc)*2*pi;

%sampling rate
fs    = 100;

%time span
time_span = [0:1/fs:75];

% initial condition
z0       = 0;    % z0 = Disp. vertikal
chi0     = 0;    % chi0= Disp. Pitch
zdot_0   = 0;    % zdot_0 = vel. vertikal
chidot_0 = 0;    % chidot_0 = vel. Pitch
%initial state vector
IC = [z0 chi0 zdot_0 chidot_0];

%% menggunakan ode45 dengan inital kondisi state vector
[t, state_vector] = ode45(@fungsi_eksitasi_harmonic,time_span,IC,'options',m,J,C1,C2,K1,K2,K3,K4,L1,L2,omega,amp,delta);
%time_span = t';
%bounce motion
z_t = state_vector(:,1);
%pitch motion
chi_t = state_vector(:,2);
%bounce velocity
v_t = state_vector(:,3);
%pitch velocity
vchi_t = state_vector(:,4);

figure(1)
subplot(1,2,1)
plot(t,z_t)
xlabel('Time[s]')
ylabel('Bounce[m]')

subplot(1,2,2)
plot(t,chi_t)
xlabel('Time[s]')
ylabel('Pitch[radians]')
