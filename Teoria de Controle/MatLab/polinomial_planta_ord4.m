clear all
clc
E = [0      0       0       0       25.72   0       0       0; 
    198.4   0       0       0       1.656   25.72   0       0;
    73.88   198.4   0       0       0       1.656   25.72   0;
    5.208   73.88   198.4   0       0       0       1.656   25.72;
    1       5.208   73.88   198.4   0       0       0       1.656;
    0       1       5.208   73.88   0       0       0       0; 
    0       0       1       5.208   0       0       0       0; 
    0       0       0       1       0       0       0       0];

A = tf([E(5,1) E(4,1) E(3,1) E(2,1) E(1,1)],1);
B = tf([E(2,5) E(1,5)],1);
G_4 = B/A

real_pole = -0.5;
zeta = 0.59115;
omega_n = -real_pole/zeta;
imag_pole = omega_n*(zeta^2-1)^(1/2);

p_1 = real_pole + imag_pole;
p_2 = real_pole - imag_pole;
p_3 = real(p_1)*10;
p_4 = p_3*10;
p_5 = p_4*10;
p_6 = p_5*10;
p_7 = p_6*10;

D = conv(conv(conv(conv(conv(conv([1 -p_1],[1 -p_2]),[1 -p_3]),[1 -p_4]),[1 -p_5]),[1 -p_6]),[1 -p_7]);
D_tf = tf(1,D);
M = E\D';

alpha = tf([M(4) M(3) M(2) M(1)],1);
beta = tf([M(8) M(7) M(6) M(5)],1);
C = beta/alpha

Gmf1 = (alpha*B)/(alpha*A + beta*B)

Gmf2 = (beta*B)/(alpha*A + beta*B)
y = pole(Gmf2);
x = zero(Gmf2);
k_P2 = (4.539*10^7)/(3.2*10^(-14)); %inverso de P2 com s->0
P2 = k_P2*tf(conv(conv(conv(conv([1 -y(3)],[1 -y(4)]),[1 -y(5)]),[1 -y(6)]),[1 -y(7)]),conv(conv(conv(conv([1 10],[1 -x(1)]),[1 -x(2)]),[1 100]),[1 1000]));

step(P2*Gmf2)
grid on
stepinfo(P2*Gmf2)