%% Limpieza
clear
close all
clc

%% Datos del problema

V=25; %[dm^3]
vo=10; %[dm^3/min]
Cao=1.5; %[mol/dm^3]
Cbo=2.0; %[mol/dm^3]
T=V/vo; %[min]

%% Reacciones

% A+2B=C+D con:
k1=0.25; %[dm^6/mol^2 min]
% 2D+3A=C+E con:
k2=0.10; %[dm^3/mol min]
% B+2C=D+F con:
k3=5.00; %[dm^6/mol^2 min]

% Leyes de Velocidad
% r1=k1*Ca*Cb^2
% r2=k2*Ca*Cd
% r3=k3*Cb*Cc^2

% Velocidades Netas de cada compuesto
% -ra=r1+3r2
% -rb=2r1+r3
% rc=r1+r2-2r3
% rd=r1+r3-2r2

%% Arreglo

% Se utilizara un reactor TAC y luego un reactor FPI con reciclo

%% Reactor TAC

% El balance para un reactor TAC es:

% Fio-Fi+ri*V=0

% Como las reacciones ocurren en fase liquida, podemos sustituir Fi=Ci*vo

% El balance quedaria como: Cio-Ci+ri*T=0 donde T=V/vo (Tao)

% Con la forma de las velocidades de reaccion, se obtiene un sistema de
% ecuaciones no lineales con los balances de A,B,C y D.

% Para resolverlo se utilizara la funcion fsolve

C_semilla=[1,1,1,1];
options=optimoptions('fsolve','Display','off');
C_TAC=fsolve(@Balance_TAC,C_semilla,options);

% Calculamos las concentraciones de E y F

Ce=T*k2*C_TAC(1)*C_TAC(4);
Cf=T*k3*C_TAC(2)*(C_TAC(3)^2);
C1=[C_TAC,Ce,Cf];

% Calculamos los flujos

F1=C1*vo;

%% Reactor FPI

% El balance para un reactor FPI es:

% dFi/dV=ri, Sustituyendo Fi=Ci*vo se obtiene: dCi/dV=ri/vo

% Se tiene que resolver un sistema de ecuaciones diferenciales, para esto
% se implementa la funcion ode45

C_inic=C1;
InterVol=[0,V];
options=odeset('RelTol',1e-4);
[Vol,C_step]=ode45(@Balance_FPI,InterVol,C_inic,options);

% Grafica del comportamiento en el FPI

plot(Vol,C_step)

% Concentraciones de salida del FPI 

C3_SR=C_step(end,:);
F3_SR=C3_SR*vo;

%% Selectividad y Rendimiento

% Selectividad D/E

Sde_SR=F3_SR(4)/F3_SR(5);

% Rendimiento

Rda_SR=F3_SR(4)/(Cao*vo-F3_SR(1)); % de D con respecto a A
Rdb_SR=F3_SR(4)/(Cbo*vo-F3_SR(2)); % de D con respecto a B

% Conversion

Xa_SR=(Cao*vo-F3_SR(1))/(Cao*vo); % Conversion de A
Xb_SR=(Cbo*vo-F3_SR(2))/(Cbo*vo); % Conversion de B

%% Reactor FPI con Reciclo

R=0.5; % Relacion de reflujo

v1=vo*(1+R); % Nuevo flujo volumetrico
F1_R=F1+F3_SR*R;
C1_R=F1_R/v1; 
Tol=1e-4;
k=0;

C_iter(1,:)=C3_SR;
F_iter(1,:)=F3_SR;
iter(1,:)=0;
error_iter(1,:)=0;

error=10;

for k=1:50
    
    [Vol_R,C_step_R]=ode45(@Balance_FPI_R,InterVol,C1_R,options);
    C3_R=C_step_R(end,:);
    F3_R=C3_R*vo;
    error=max(F1_R-R*F3_R-F1);
    F1_R=F1+R*F3_R;
    C1_R=F1_R/v1;
    
    
    C_iter(k+1,:)=C3_R;
    F_iter(k+1,:)=F3_R;
    iter(k+1)=k;
    error_iter(k+1)=error;
    if error<=Tol
        break
    end
end

%% Selectividad y Rendimiento con Reciclo

% Selectividad D/E

Sde_R=F3_R(4)/F3_R(5);

% Rendimiento

Rda_R=F3_R(4)/(Cao*vo-F3_R(1)); % de D con respecto a A
Rdb_R=F3_R(4)/(Cbo*vo-F3_R(2)); % de D con respecto a B

% Conversion

Xa_R=(Cao*vo-F3_R(1))/(Cao*vo); % Conversion de A
Xb_R=(Cbo*vo-F3_R(2))/(Cbo*vo); % Conversion de B

%% Selectividades instantaneas

% Velocidades

r1=@(C)k1*C(1)*(C(2)^2);
r2=@(C)k2*C(1)*C(4);
r3=@(C)k3*C(2)*(C(3)^2);

rd=@(C)r1(C)+r3(C)-2*r2(C);
re=@(C)r2(C);

% Selectividad instantanea a la salida del TAC

Sde_ITAC=rd(C1)/re(C1);

% Selectividad instantanea a la salida del FPI sin Reciclo

Sde_IFPI=rd(C3_SR)/re(C3_SR);

% Selectividad instantanea a la entrada FPI con Reciclo

Sde_IER=rd(C1_R)/re(C1_R);

% Selectividad instantanea a la salida FPI con Reciclo

Sde_IFPIR=rd(C3_R)/re(C3_R);