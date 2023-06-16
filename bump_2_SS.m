clear all;
clc;
K1    =580000;
K2    =580000; 
K3    =580000; 
K4    =580000; %Newton
m     =29892; 
J     =29892; %Newton
C1    =34000;
C2    =34000; %Newton

%Panjang titik pusat ke titik pusat roda
L1    =1.5; 
L2    =1.5; %Meter

r=4;
v=15*10^3; %kecepatan awal

%%Matrix MJ C K F
M       =[m 0 ; 0 J]; %#ok<GPFST> %%MASS MATRIX
C       =[C1+C2 -C1*L1+C2*L2; -C1*L1+C2*L2 C1*L1^2+C2*L2^2]; %%STIFFNES MATRIX
K       =[K1+K2+K3+K4 -K1*L1+K2*L2-K3*L1+K4*L2; -K1*L1+K2*L2-K3*L1+K4*L2 K1*L1^2+K2*L2^2+K3*L1^2+K4*L2^2]; %%DAMPED MATRIX
kz      =[K1+K3 K2+K4 ; -K1*L1-K3*L1 K2*L2+K4*L2]; %Eksitasi
cz_dot  =[C1 C2; -C1*L1 C2*L2]; %Eksitasi
F       =[kz + cz_dot]; %#ok<GPFST>

A00     =zeros(2);
A11     =eye(2);


%%eksitasi bump
sim = 10;
c = 0
a1 = 0.5;
a2 = 0.3;
for i= 0:0.01:sim;
    c=c+1;
    if i<0
        ud(c)=a1*(1-cos(4*pi*i));
    elseif i>2 && i<2.5
        ud(c)=a2*(1-cos(4*pi*i));
    else
        ud(c)=0;
    end
end
i = 0:0.01:sim;
plot(i,ud,'Linewidth',2);

A  = [A00 A11; -inv(M)*K -inv(M)*C]; %#ok<GPFST>
B  = [A00 ; inv(M)*F];
C  = [1 1 0 0];
D  = [0 0];

SYS = ss(A,B,C,D);
%%simulasi
tsim = 0:0.01:sim;
uc = zeros(size(tsim));
u  = [uc; ud];
x0 = [0 0 0 0];
y  = lsim(SYS,u,tsim, x0);

i = 0:0.01:sim;
plot(i,y,'b',i,ud,'k','Linewidth',2);
