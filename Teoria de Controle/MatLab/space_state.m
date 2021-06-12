clear all
close all
clc
offset = 2;
wt = 4;

% ESPECIFICAÇÔES DE DESEMPENHO %
OS = 10;
ts = 4;
zeta = -log(OS/100)/(sqrt(pi^2 + (log(OS/100) )^2));
wn = 4/(zeta*ts);
Dom_poles = [-zeta*wn + wn*sqrt(zeta^2-1) -zeta*wn - wn*sqrt(zeta^2-1)];
P = [Dom_poles real(Dom_poles(1)*10)]
%n+1 pólos graças à presença da ação integral na matriz Aa

% FUNÇÃO DE TRANSFERÊNCIA E ESPAÇO DE ESTADOS %
load('C:\Users\Maurício\Google Drive\CEFET\7° Período\Teoria de Controle (60h)\Planta\validacaoOK_2.mat')

[A,B,C,D] = tf2ss(G.Numerator{1},G.Denominator{1});

% REALIMENTAÇÃO DE ESTADOS COM AÇÃO INTEGRAL %
%Confere se o sistema é controlável
if(size(A)==rank(ctrb(A,B)))
    'Controlável'
end

%%
Aa = [A zeros(size(A,1),1);-C 0]; %A aumentada pela ação integral
Ba = [B; 0]; %B aumentada pela ação integral
K_ = ones(1, size(Aa,1));
Fctr = [real(P(1)) imag(P(1)) 0; imag(P(2)) real(P(2)) 0; 0 0 real(P(3))]
%Fctr = [real(P(1)) imag(P(1)) 0 0; imag(P(2)) real(P(2)) 0 0; 0 0 zero(G) 0; 0 0 0 real(P(3))]
%Cancelando o zero pela alocação de polos

Tctr = lyap(Aa, -Fctr, -Ba*K_);
Kt = K_/Tctr;
K = Kt(1:size(A,1))
Ka = -Kt(end)   

%%
% OBSERVADOR ULTRAVELOZ %
%Confere se o sistema é observável
if(size(A)==rank(obsv(A,C)))
    'Observável'
end

L = ones(size(A,1),1);
%Fobs = diag([real(P(1))*5 real(P(1))*5-1 real(P(1))*5-2])
%Fobs = diag([real(P(1))*10 real(P(1))*10-1 real(P(1))*10-2])
Fobs = diag([real(P(1))*15 real(P(1))*15-1])

if(size(A)==rank(ctrb(Fobs,L)))
    'Controlável'
end

Tobs = lyap(-Fobs, A, -L*C);
% Acompanha a din   âmica, independentemente do ruído branco.
% Contudo, o sinal de controle fica mais suscetível à saturação %