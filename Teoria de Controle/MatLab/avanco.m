clear all;
close all;
clc

load('C:\Users\Maurício\Google Drive\CEFET\7° Período\Teoria de Controle (60h)\Planta\validacaoOK_2.mat')

syms s
[Num,Den] = tfdata(G);
G_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);

OS = 10; %10
ts = 4; %4
zeta = -log(OS/100)/(sqrt(pi^2 + (log(OS/100) )^2));
wn = 4/(zeta*ts);
P = [-zeta*wn + wn*sqrt(zeta^2-1) -zeta*wn - wn*sqrt(zeta^2-1)];
poles = pole(G);

figure
plot(zpk(G).P{1}, [0 0], 'rx');
grid
hold on;
plot(real(P), imag(P), 'bs');

alpha_1 = 180 - atand(imag(P(1))/-real(P(1)));
alpha_2 = atand(imag(P(1))/(-poles(2)+real(P(1))));
phi = 180 - (alpha_1 + alpha_2);

beta = 90 + atand(-real(P(1))/imag(P(1)));
bissetriz = beta/2;

phi_p = bissetriz - phi/2 - atand(-real(P(1))/imag(P(1)));
phi_z = phi + phi_p;

p = -real(P(1)) + tand(phi_p)*imag(P(1));
z = -real(P(1)) + tand(phi_z)*imag(P(1));

C = (s + z)/(s + p);
Kc = 1/abs(subs(G_syms*C,'s',P(1)));
C = tf(double(vpa(Kc))*[1 z],[1 p])

plot(-z, 0, 'go')
plot(-p, 0, 'gx')

Gmf = feedback(C*G,1);
Filtro = tf([1 -zpk(Gmf).P{1, 1}(1)],[1 -zero(C)]);
Filtro = Filtro/dcgain(Filtro);

stepinfo(Filtro*Gmf)
