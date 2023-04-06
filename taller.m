%% Script Principal

% Reactores con reciclo.
%Esmeralda López 12-10576
%Kevin Mejias 11-10600
%Gerardo Hernandez 12-10382


clc; clear all; %Limpiar pantalla
format short g;  %Formato

%% Datos
Tau=2.5;  %[min] tiempo de residencia
vo=10;    %[L/min] caudal de entrada al TAC
V=25;     %[dm^3] volumen del TAC y del FPI
 
k1=0.25;  %[dm^6/mol^2min] constante de velocidad 1 
k2=0.10;  %[dm^3/molmin] constante de velocidad 2
k3=5;     %[dm^6/mol^2min] constante de velocidad 3

CAO=1.5;  %[mol/L] concentracion de A que entra al TAC
CBO=2;    %[mol/L] concentracion de B que entra al TAC

%% Estimado inicial de concentraciones
CO=[1;1;1;1;1;1];
%% Creacion funcion anonima
fun=@(C)TAC(C,CAO,CBO,k1,k2,k3,Tau);
%% Resolucion
[sol,fval,bandera,iterTAC] = fsolve(fun,CO);
CA=sol(1);
CB=sol(2);
CC=sol(3);
CD=sol(4);
CE=sol(5);
CF=sol(6);
iterTAC;

fprintf('Los flujos a la salida del TAC son \b\n')
FATAC=CA*vo
FBTAC=CB*vo
FCTAC=CC*vo
FDTAC=CD*vo
FETAC=CE*vo
FFTAC=CF*vo

%% FPI sin reciclo 
FO=[FATAC;FBTAC;FCTAC;FDTAC;FETAC;FFTAC];
V=linspace(0,V);
fun=@(F)FPI(F,V,vo);
[F,istate,msg]=lsode(fun,FO,V,vo);

fprintf('Los flujos a la salida del FPI sin reciclo son \b\n')
FAFPI=F(end,1)
FBFPI=F(end,2)
FCFPI=F(end,3)
FDFPI=F(end,4)
FEFPI=F(end,5)
FFFPI=F(end,6)

fpi=[FAFPI;FBFPI;FCFPI;FDFPI;FEFPI;FFFPI];
S=selectividad(k1,k2,k3,fpi,vo);
fprintf('La selectividad de D con respecto a E es %f \b\n',S)

%% Resolviendo con un reciclo R de 50% en 3 iteraciones
R=0.5;
v=vo*(1+R);
iter=0;
FO=fpi;

for i=1:3
fun=@(F)FPI(F,V,v);
[F,istate,msg]=lsode(fun,FO,V,v);
FAFPI=F(end,1);
FBFPI=F(end,2);
FCFPI=F(end,3);
FDFPI=F(end,4);
FEFPI=F(end,5);
FFFPI=F(end,6);

FA=FATAC+(R*FAFPI);
FB=FBTAC+(R*FBFPI);
FC=FCTAC+(R*FCFPI);
FD=FDTAC+(R*FDFPI);
FE=FETAC+(R*FEFPI);
FF=FFTAC+(R*FFFPI);
FO=[FA;FB;FC;FD;FE;FF];
iter=iter+1;
end

fprintf('Los flujos a la salida del FPI con reciclo despues de 3 iteraciones son \b\n')
FAFPI=F(end,1)
FBFPI=F(end,2)
FCFPI=F(end,3)
FDFPI=F(end,4)
FEFPI=F(end,5)
FFFPI=F(end,6)

fprintf('Los flujos a la entrada del FPI con reciclo despues de 3 iteraciones son \b\n')
FA
FB
FC
FD
FE
FF

fpi3it=[FAFPI;FBFPI;FCFPI;FDFPI;FEFPI;FFFPI];
S=selectividad(k1,k2,k3,fpi3it,vo);
fprintf('La selectividad de D con respecto a E despues de 3 iteraciones es %f \b\n',S)

%% Resolviendo con reciclo de 50% para una tolerancia
FO=fpi;
tol=0.0001;
n=1000;
iter=0;
Fs=fpi;

for i=1:n
fun=@(F)FPI(F,V,v);
[F,istate,msg]=lsode(fun,FO,V,v);
FAFPI=F(end,1);
FBFPI=F(end,2);
FCFPI=F(end,3);
FDFPI=F(end,4);
FEFPI=F(end,5);
FFFPI=F(end,6);

FA=FATAC+R*FAFPI;
FB=FBTAC+R*FBFPI;
FC=FCTAC+R*FCFPI;
FD=FDTAC+R*FDFPI;
FE=FETAC+R*FEFPI;
FF=FFTAC+R*FFFPI;
aux=[FAFPI;FBFPI;FCFPI;FDFPI;FEFPI;FFFPI];
FO=[FA;FB;FC;FD;FE;FF];
iter=iter+1;

if any((aux-Fs)>tol)
Fs=aux;
elseif all((aux-Fs)<tol)
fprintf('Se logro satisfacer una tolerancia de 0.0001 despues de %d iteraciones \b\n',iter)
fprintf('Los flujos a la salida del FPI con reciclo satisfaciendo la tolerancia son \b\n')
FAFPI=F(end,1)
FBFPI=F(end,2)
FCFPI=F(end,3)
FDFPI=F(end,4)
FEFPI=F(end,5)
FFFPI=F(end,6)
fprintf('Los flujos a la entrada del FPI con reciclo satisfaciendo la tolerancia son \b\n')
FA
FB
FC
FD
FE
FF
break
end
end

fpifi=[FAFPI;FBFPI;FCFPI;FDFPI;FEFPI;FFFPI];
S=selectividad(k1,k2,k3,fpifi,vo);
fprintf('La selectividad de D con respecto a E satisfaciendo la tolerancia es %f \b\n',S)