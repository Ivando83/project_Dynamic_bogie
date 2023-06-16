function [dwdt] = fungsi_eksitasi_harmonik(t,w,m,J,C1,C2,K1,K2,K3,K4,L1,L2,omega,amp,delta)
w1 = w(1);%bounce motion
w2 = w(2);%pitch motion
w3 = w(3);%bounce velocity
w4 = w(4);%pitch velocity

%bounce velocity
dw1dt = w3;
%pitch velocity
dw2dt = w4;
%bounce acceleration
dw3dt = -((K1+K2+K3+K4)/m)*w1 - ((-K1*L1-K3*L1+K2*L2+K4*L2)/m)*w2 - ((C1+C2)/m)*w3 - ((-C1*L1+C2*L2)/m)*w4 + ((K1*amp*cos(omega*t))/m) + ...
    ((K3*amp*cos(omega*t))/m) - ((C1*amp*omega*sin(omega*t))/m) + ((K2*amp*cos(omega*t-delta))/m) + ((K4*amp*cos(omega*t-delta))/m) - ((C2*amp*omega*sin(omega*t-delta))/m);
%pitch acceleration
dw4dt = -((-K1*L1-K3*L1+K2*L1+K4*L2)/J)*w1 - ((K1*(L1)^2+K2*(L2)^2+K3*(L1)^2+K4*(L2)^2)/J)*w2 - ((-C1*L1+C2*L2)/J)*w3 - ((C1*(L1)^2+C2*(L2)^2)/J)*w4 + ...
    ((-K1*L1*amp*cos(omega*t))/J) - ((K3*L1*amp*cos(omega*t))/J) + ((C1*L1*amp*sin(omega*t))/J) + ((K2*L2*amp*cos(omega*t-delta))/J) + ((K4*L2*amp*cos(omega*t-delta))/J)...
    - ((C2*L2*omega*amp*sin(omega*t-delta))/J);

%first derivative of the state matrix
dwdt = [dw1dt;dw2dt;dw3dt;dw4dt];



