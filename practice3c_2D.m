
%%%%%%%%%%%%%%%%%%%%% Program Persamaan Gerak dengan Eksitasi Eksitasi Bump Modifikasi %%%%%%%%%%%%%%%%%% Persamaan Gerak dengan eksitasi Bump
K1    =580000;     %koefisien spring depan [N/m]
K2    =580000;     %koefisien spring belakang [N/m]
K3    =580000;     %koefisien spring depan [N/m]
K4    =580000;     %koefisien spring belakang [N/m]
m     =3048.1357;  %Massa Bogie [Kg]
J     =280.4219;   %Inertia Bogie [Kg*m^2]
C1    =34000;      %koefisien Damping depan [Ns/m]
C2    =34000;      %koefisien Damping belakang [Ns/m]

%Panjang titik pusat ke titik pusat roda
L1    = 1.5; %Meter
L2    = 1.5; %Meter
v     = 16;  %kecepatan awal m/s
r     = 0.008; %jarijari m

%%Matrix M C K F
% Mass Matrix
M       =[m 0 ; 0 J]; %#ok<GPFST> %%MASS MATRIX
% Damping Matrix
C       =[C1+C2 -C1*L1+C2*L2; -C1*L1+C2*L2 C1*L1^2+C2*L2^2]; 
% Stiffnes Matrix
K       =[K1+K2+K3+K4 -K1*L1+K2*L2-K3*L1+K4*L2;...
         -K1*L1+K2*L2-K3*L1+K4*L2 K1*L1^2+K2*L2^2+K3*L1^2+K4*L2^2]; 
     

% Force Excitatin for Stiffnes and damping
kz      = [K1+K3 K2+K4 ; -K1*L1-K3*L1 K2*L2+K4*L2]; %Eksitasi Stiffnes
cz      = [C1 C2; -C1*L1 C2*L2]                     %Eksitasi demping

%%Force F(t) Input Bump
f       = @(time)  sqrt(r^2-(time-r)^2);    
fdot    = @(time)  sqrt(r^2-(v.*time-r)^2); 

A0 = zeros(2);
I  = eye(2);
%% State Variable
A0 = zeros(2);
I  = eye(2);
A = [A0 I; -inv(M)*K -inv(M)*C]; %Matrix A
B = [A0; inv(M)]; %Matrix B
F = @(TIME)[[kz]*[ f(TIME); 0]+[cz]*[fdot(TIME); 0]]; %Matrix Eksitasi

%%State Space sdot = A*s+B*F
sdot = @(t,s) A*s(1:4) + B*F(t);

% initial condition
z0       = 0;    % z0 = Disp. vertikal
chi0     = 0;    % chi0= Disp. Pitch
zdot_0   = 0;    % zdot_0 = vel. vertikal
chidot_0 = 0;    % chidot_0 = vel. Pitch
IC = [z0 chi0 zdot_0 chidot_0];

%Time span 
t0=0; tf=2;
tspan = [t0:0.01:tf]; %Rentang Waktu

%%Numerik
[time state_values]=ode45(sdot, tspan, IC);
z       = state_values(:,1); %disp bounch
chi     = state_values(:,2); %disp pitch
zdot    = state_values(:,3); %vel bounch
chidot  = state_values(:,4); %vel pitch

%%plot
figure(1),clf
plot(time,z), xlabel('Time (s)'), ylabel('Displacment (m)')
title ('Respons bounce displacment vs Time')

figure(2),clf
plot(time,chi), xlabel('Time (s)'), ylabel('Displacment [radians]')
title('Respons Pitch Roda Belakang vs Time')

figure(3),clf
plot(time,zdot), xlabel('Time (s)'), ylabel('Vertical Velocity')
title('Kecepatan bounce vs Time')

figure(4),clf
plot(time,chidot), xlabel('Time (s)'), ylabel('Vertical Velocity')
title('Kecepatan Pitch vs Time')





