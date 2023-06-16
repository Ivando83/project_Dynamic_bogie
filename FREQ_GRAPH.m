%plotting the frequency repsonse graph
clear all %#ok<CLALL>
close all
clc

K1    =580000; %Newton
K2    =580000; %Newton
K3    =580000; %Newton
K4    =580000; %Newton
m     =29892;  %Newton
J     =2750;  %Newton
C1    =34000;  %Newton
C2    =34000;  %Newton

%Panjang titik pusat ke titik pusat roda
L1    =1.5;    %Meter
L2    =1.5;    %Metere
v     =16;     %Kecepatan awal m/s

panjang_eksitasi = 0.55;  %Panjang antara bantalan (m)
omeg  = (v*2*pi)/panjang_eksitasi;     %Frekuensi (Hz)
delta = ((L1+L2)*2*pi)/panjang_eksitasi; %Sudut fase
amp   = 8*10^-3;                         %Amplitido

omega = 0:0.01:omeg;
tspan = 0:0.1:8;
%z = zeros(2,size(omega,2));

for kk = 1:length(omega)
    %%Matrix M C K F
% Mass Matrix
M       =[m 0 ; 0 J]; %#ok<GPFST> %%MASS MATRIX
% Damping Matrix
C       =[C1+C2 -C1*L1+C2*L2; -C1*L1+C2*L2 C1*L1^2+C2*L2^2]; 
% Stiffnes Matrix
K       =[K1+K2+K3+K4 -K1*L1+K2*L2-K3*L1+K4*L2; ...
    -K1*L1+K2*L2-K3*L1+K4*L2 K1*L1^2+K2*L2^2+K3*L1^2+K4*L2^2]; 

% Eksitasi Matriks 
% h_c (cosinus)
h_c = [(C2*amp*omega(kk)*sin(delta) + K1*amp + K2*amp*cos(delta)...
    + K3*amp + K4*amp*cos(delta))*cos(omega(kk)*tspan);...
    (C2*L2*amp*omega(kk)*sin(delta) - K1*L1*omega(kk) + K2*L2*omega(kk)*cos(delta)...
    - K3*L1*omega(kk) + K4*L2*omega(kk)*cos(delta))*cos(omega(kk)*tspan)];

% h_s (sinus)
h_s = [(-C1*amp*omega(kk) - C2*amp*omega(kk)*cos(delta) + K2*amp*sin(delta) ...
    + K4*amp*sin(delta))*sin(omega(kk)*tspan);...
    (C1*L1*amp*omega(kk) - C2*L2*amp*omega(kk)*cos(delta) + K2*L2*amp*sin(delta) ...
    + K4*L2*amp*sin(delta))*sin(omega(kk)*tspan)];

h_star = 0.5*(h_c - 1i*h_s);

F_star = inv(-(omega(kk)^2)*M + 1i*omega(kk)*C + K); %2x2 matrix

x_star = F_star*h_star; %#ok<MINV>


%bounc cos and sin particular 
z_c = real(2*x_star(1));
z_s = -imag(2*x_star(1));

%z(:,kk) = sqrt(z_c.^2 + z_s.^2);

%pitch
chi_c = real(2*x_star(2));
chi_s = -imag(2*x_star(2));

%amplitude of the motion for the particular omega
%bounce
z_amp(kk) = sqrt((z_c)^2 + (z_s)^2);
%pitch
chi_amp(kk) = sqrt((chi_c)^2 + (chi_s)^2);
end

figure(1)
subplot(1,2,1)
plot(omega,z_amp)
title('Bounce Frequency Response')
xlabel('Frekuensi (Hz)')
ylabel('Bounce Amplitude (m)')

subplot(1,2,2)
plot(omega,chi_amp)
title('Pitch Frequency Response')
xlabel('Frekuensi (Hz)')
ylabel('Pitch Amplitude ')


%z_hat = z(1,:);
%chi_hat = z(2,:);


