%% RQ_1110600_1210576_1210382
% Reactores con reciclo.
clc; clear all; %Limpiar pantalla
format short g;  %Formato
%% Inicialmente, definimos las variables a introducir:
% Datos generales del reactor:
k1=0.25;
k2=0.1;
k3=5.00;
V=25;
v0=10;
CA0=1.5;
CB0=2.0;
Tao=V/v0;

%% Ecuaciones a resolver

% ra1=-k1*CA*CB^2
% rb1=-2*k1*CB^2
% rc1=k1*CA*CB^2
% rd1=k1*CA*CB^2
% ra2=-3*k2*CA*CD
% rc2=k2*CA*CD
% rd2=-2*k2*CA*CD
% re2=k2*CA*CD
% rb3=-k3*CB*CC^2
% rc3=-2*k3*CB*CC^2
% rd3=k3*CB*CC^2
%rf3=k3*CB*CC^2

% ra=ra1+ra2
% rb=rb1+rb3
% rc=rc1+rc2+rc3
% rd=rd1+rd2+rd3
% re=re2
% rf=rf3

%Sustituyendo en las velocidades totales

% ra=-k1*CA*CB^2-3*k2*CA*CD
% rb=-2*k1*CB^2-k3*CB*CC^2
% rc=k1*CA*CB^2+k2*CA*CD-2*k3*CB*CC^2
% rd=k1*CA*CB^2-2*k2*CA*CD+k3*CB*CC^2
% re=k2*CA*CD
% rf=k3*CB*CC^2

% %% Balances
% CA0-CA+ra*Tao=0
% CB0-CB+rb*Tao=0
% -CC+rc*Tao=0
% -CD+rd*Tao=0
% -CE+re*Tao=0
% -CF+rf*Tao=0

% Sustituimos para obtener 6 Ec con 6 incognitas

% CA0-CA+(-k1*CA*CB^2-3*k2*CA*CD)*Tao=0
% CB0-CB+(-2*k1*CB^2-k3*CB*CC^2)*Tao=0
% -CC+(k1*CA*CB^2+k2*CA*CD-2*k3*CB*CC^2)*Tao=0
% -CD+(k1*CA*CB^2-2*k2*CA*CD+k3*CB*CC^2)*Tao=0
% -CE+(k2*CA*CD)*Tao=0
% -CF+(k3*CB*CC^2)*Tao=0

%Construimos la matriz de las ecuaciones a resolver

F=@(C) [CA0-C(1)+(-k1*C(1)*C(2)^2-3*k2*C(1)*C(4))*Tao;
CB0-C(2)+(-2*k1*C(2)^2-k3*C(2)*C(3)^2)*Tao;
-C(3)+(k1*C(1)*C(2)^2+k2*C(1)*C(4)-2*k3*C(2)*C(3)^2)*Tao;
-C(4)+(k1*C(1)*C(2)^2-2*k2*C(1)*C(4)+k3*C(2)*C(3)^2)*Tao;
-C(5)+(k2*C(1)*C(4))*Tao;
-C(6)+(k3*C(2)*C(3)^2)*Tao];
C0=[8, 8, 1, 4, 1, 2];
C0=fsolve(F,C0);
F1=C0*v0;
F1=F1';

%%Ahora realizamos los balances para el FPI
%dFA/dV=ra
%dFB/dV=rb
%dFC/dV=rc
%dFD/dV=rd
%dFE/dV=re
%dFF/dV=rf
%Sustituyendo los valores de las velocidades
%dFA/dV=-k1*CA*CB^2-3*k2*CA*CD
%dFB/dV=-2*k1*CB^2-k3*CB*CC^2
%dFC/dV=k1*CA*CB^2+k2*CA*CD-2*k3*CB*CC^2
%dFD/dV=k1*CA*CB^2-2*k2*CA*CD+k3*CB*CC^2
%dFE/dV=k2*CA*CD
%dFF/dV=k3*CB*CC^2
F=[4.0227 4.7927 1.0990 5.8995 1.6695]
cond=[FA(0)==F1(1), FB(0)==F1(2), FC(0)==F1(3), FD(0)==F1(4), FE(0)==F1(5), FF(0)==F1(6)];
F2=dsolve('DFA=-k1*CA*CB^2-3*k2*CA*CD','DFB=-2*k1*CB^2-k3*CB*CC^2','DFC=k1*CA*CB^2+k2*CA*CD-2*k3*CB*CC^2','DFD=k1*CA*CB^2-2*k2*CA*CD+k3*CB*CC^2','DFE=k2*CA*CD','DFF=k3*CB*CC^2','V',cond) 
FA=F2.FA;
FB=F2.FB;
FC=F2.FC;
FD=F2.FD;
FE=F2.FE;
FE=F2.FE;


