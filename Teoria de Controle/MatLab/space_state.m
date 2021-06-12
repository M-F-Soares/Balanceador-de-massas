clear all
close all
clc
offset = 2;
wt = 4;

% ESPECIFICA��ES DE DESEMPENHO %
OS = 10;
ts = 4;
zeta = -log(OS/100)/(sqrt(pi^2 + (log(OS/100) )^2));
wn = 4/(zeta*ts);
Dom_poles = [-zeta*wn + wn*sqrt(zeta^2-1) -zeta*wn - wn*sqrt(zeta^2-1)];
P = [Dom_poles real(Dom_poles(1)*10)]
%n+1 p�los gra�as � presen�a da a��o integral na matriz Aa

% FUN��O DE TRANSFER�NCIA E ESPA�O DE ESTADOS %
load('C:\Users\Maur�cio\Google Drive\CEFET\7� Per�odo\Teoria de Controle (60h)\Planta\validacaoOK_2.mat')

[A,B,C,D] = tf2ss(G.Numerator{1},G.Denominator{1});

% REALIMENTA��O DE ESTADOS COM A��O INTEGRAL %
%Confere se o sistema � control�vel
if(size(A)==rank(ctrb(A,B)))
    'Control�vel'
end

%%
Aa = [A zeros(size(A,1),1);-C 0]; %A aumentada pela a��o integral
Ba = [B; 0]; %B aumentada pela a��o integral
K_ = ones(1, size(Aa,1));
Fctr = [real(P(1)) imag(P(1)) 0; imag(P(2)) real(P(2)) 0; 0 0 real(P(3))]
%Fctr = [real(P(1)) imag(P(1)) 0 0; imag(P(2)) real(P(2)) 0 0; 0 0 zero(G) 0; 0 0 0 real(P(3))]
%Cancelando o zero pela aloca��o de polos

Tctr = lyap(Aa, -Fctr, -Ba*K_);
Kt = K_/Tctr;
K = Kt(1:size(A,1))
Ka = -Kt(end)   

%%
% OBSERVADOR ULTRAVELOZ %
%Confere se o sistema � observ�vel
if(size(A)==rank(obsv(A,C)))
    'Observ�vel'
end

L = ones(size(A,1),1);
%Fobs = diag([real(P(1))*5 real(P(1))*5-1 real(P(1))*5-2])
%Fobs = diag([real(P(1))*10 real(P(1))*10-1 real(P(1))*10-2])
Fobs = diag([real(P(1))*15 real(P(1))*15-1])

if(size(A)==rank(ctrb(Fobs,L)))
    'Control�vel'
end

Tobs = lyap(-Fobs, A, -L*C);
% Acompanha a din   �mica, independentemente do ru�do branco.
% Contudo, o sinal de controle fica mais suscet�vel � satura��o %