close all;
clear all;
clc;

E = [0      0       0       3.502   0       0; 
     45.804 0       0       0       3.502   0;
     3.5045 45.804  0       0       0       3.502;
     1      3.5045  45.804  0       0       0;
     0      1       3.5045  0       0       0;
     0      0       1       0       0       0];
 
A = tf([E(4,1) E(3,1) E(2,1) E(1,1)],1);
B = tf(E(1,4),1);
G = B/A

real_pole = -0.25;
zeta = 0.357857130503;
omega_n = -real_pole/zeta;
imag_pole = omega_n*(zeta^2-1)^(1/2);

p_1 = real_pole + imag_pole;
p_2 = real_pole - imag_pole;
p_3 = real(p_1)*10;
p_4 = p_3*10;
p_5 = p_4*10;

D = conv(conv(conv(conv([1 -p_1],[1 -p_2]),[1 -p_3]),[1 -p_4]),[1 -p_5]);
D_tf = tf(1,D);
D = flip(D');
M = E\D;

alpha = tf([M(3) M(2) M(1)],1);
beta = tf([M(6) M(5) M(4)],1);
C = beta/alpha

Gmf1 = (alpha*B)/(alpha*A + beta*B)%paralelo
y = pole(Gmf1);
x = zero(Gmf1);
k_P1(1) = 1/(Gmf1.Numerator{1}(end)/Gmf1.Denominator{1}(end)); %inverso de P1 com s->0
P1 = tf(conv(conv([1 -y(1)],[1 -y(2)]),[1 -y(3)]),conv(conv([1 -x(1)],[1 -x(2)]),[1 2502]));
k_P1(2) = 1/(P1.Numerator{1}(end)/P1.Denominator{1}(end));
P1 = k_P1(1)*P1*k_P1(2);

figure, step(P1*Gmf1)
grid on
stepinfo(P1*Gmf1)


Gmf2 = (beta*B)/(alpha*A + beta*B)%série
y = pole(Gmf2);
x = zero(Gmf2);
k_P2(1) = 1/(Gmf2.Numerator{1}(end)/Gmf2.Denominator{1}(end)); %inverso de P1 com s->0
P2 = tf(conv(conv([1 -y(1)],[1 -y(2)]),[1 -y(3)]),conv(conv([1 -x(1)],[1 -x(2)]),[1 40]));
k_P2(2) = 1/(P2.Numerator{1}(end)/P2.Denominator{1}(end));
P2 = k_P2(1)*P2*k_P2(2);

figure, step(P2*Gmf2)
grid on
stepinfo(P2*Gmf2)

poles_G = pole(G);
poles_Gmf_x = [real_pole real_pole];
poles_Gmf_y = [imag_pole/j imag_pole*j];

figure 
rlocus(G)
hold on
grid

if (abs(real_pole) <= abs(real(poles_G(2))))
    F1 = tf([1],[1/(-10*(real(poles_G(2)))) 1])
    F = F1;
else
    F2 = tf([1],[1/(-10*(real_pole)) 1])
    F = F2;
end

plot([poles_Gmf_x pole(F)],[poles_Gmf_y 0],'o')