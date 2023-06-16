function [dwdt] = fungsi_eksitasi_harmonic(time_span,w,m,J,C1,C2,K1,K2,K3,K4,L1,L2,omega,amp,delta)
A = [0 0 1 0; 0 0 0 1;...
    -(-K1+K2+K3+K4)/m,-(K1*L1+K2*L2-K3*L1+K4*L2)/m,-(C1+C2)/m,-(C1*L1+C2*L2)/m;...
    -(K1*L1+K2*L1-K3*L1+K4*L2)/J,-(K1*(L1)^2+K2*(L2)^2)/J,-(C1*L1+C2*L2)/J,-(C1*(L1)^2+C2*(L2)^2)/J];
% Set up force 1
F1 = K1*cos(omega*time_span)+K3*cos(omega*time_span)-C1*omega*sin(omega*time_span)+K2*cos(omega*time_span-delta)+K4*cos(omega*time_span-delta)-C2*omega*sin(omega*time_span-delta);

% Set up force 2
F2 = (-K1*cos(omega*time_span)-K3*cos(omega*time_span)+C1*omega*sin(omega*time_span))*L1+(K2*cos(omega*time_span-delta)+K4*cos(omega*time_span-delta)-C2*omega*sin(omega*time_span-delta))*L2;
% Combine force 1 and force 2 to excitation vector
B = [0; 0; (F1*amp)/m; (F2*amp)/J];

% Create state space representation
dwdt = A*w+B;