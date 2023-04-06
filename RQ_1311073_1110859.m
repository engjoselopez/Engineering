%% RQ_1311073_1110859
% Reactores con reciclo.
%Carlos Pereira   %Carnet. 13_11073
%Oscar Reyes      %Carnet. 11_10859

clc; clear all; %Limpiar pantalla
format short g;  %Formato

%% Inicialmente, definimos las variables a introducir:
global k1 k2 vo Cao Vtac Fao tao tao2

% Datos generales del reactor:
k1=0.15;
k2=0.06;
vo=5;
Cao=0.1;
Vtac=10;
Fao=Cao*vo;
tao=Vtac/vo;
R=0.5;
 tipo_RQ=1; %(Tipo TAC)
% tipo_RQ=2; %(Tipo FPI)

if tipo_RQ==1
    % Balances de Masa para un TAC (A->B->C)
    %Comp A, B y C:
    Fa1=(Fao)/(1+k1*tao);
    Fb1=(k1*Fa1*tao)/(1+k2*tao);
    Fc1=k2*Fb1*tao;   
    %Flujos sin reciclo como flujos de entrada
    Fa0=Fa1;
    Fb0=Fb1;
    Fc0=Fc1;
    tol=0.1;
    
    %TAC con reciclo
    v=vo*(1+R);
    tao2=Vtac/v;    
    %resuelvo el TAC por primera vez con reciclo
    Fa2=Fa0/(1+k1*tao2);
    Fb2=(Fb0+Fa2*(k1*tao2))/(1+k2*tao2);
    Fc2=Fc0+(k2*tao2*Fb2);
    Fi=[Fa2,Fb2,Fc2];       %Vector de flujos iniciales(con reciclo)
 
    Flujos_iter=[Fi];  
    k=1;
    while tol>0.01

        %Ecuacion de reciclo:
        Fa1=Fa0+(R*Fa2);
        Fb1=Fb0+(R*Fb2);
        Fc1=Fc0+(R*Fc2);
        
        %Nuevos Flujos de entrada con reciclo
        Fa0=Fa1;
        Fb0=Fb1;
        Fc0=Fc1;
        
        %Resolucion del TAC por segunda vez
        Fa2=Fa0/(1+k1*tao2);
        Fb2=(Fb0+Fa2*(k1*tao2))/(1+k2*tao2);
        Fc2=Fc0+(k2*tao2*Fb2);
        Fk=[Fa2,Fb2,Fc2];              %Vector de flujos iter.k(con reciclo)
        Flujos_iter=[Flujos_iter; Fk]; %Tabla de flujos para cada iteracion.    
 
        %Calculo tolerancia
        tol=max(abs((Flujos_iter(k,:)-Flujos_iter(k+1,:))));
        k=k+1;
        if k>=20
            break
        end
    end
    Convergencia_iteracion=k
    Flujos_iter
    Flujos_salida=Flujos_iter(end,:)

    hold on
    plot([1:k],Flujos_iter(:,1),'*b')
    plot([1:k],Flujos_iter(:,2),'og')
    plot([1:k],Flujos_iter(:,3),'-y')
    
elseif tipo_RQ==2
    % Balances de Masa para un FPI(A->B->C)
    %Comp A, B y C:
    %@SED_1311073_1110859 BALANCES DE A,B y C. REACTOR FPI.
    FiInicial=[Fao 0 0];
    [Volumen,Fi1]=ode45('SED_1311073_1110859',[0 10],FiInicial);
        %Flujos sin reciclo como flujos de entrada
    Fa0=Fi1(end,1);
    Fb0=Fi1(end,2);
    Fc0=Fi1(end,3);
    Fi0=[Fa0 Fb0 Fc0]
    tol=0.1;
    
    %FPI con reciclo
    v=vo*(1+R);
    tao2=Vtac/v;    
    %resuelvo el FPI por primera vez con reciclo
    [Volumen,Fi2]=ode45('SED_1311073_1110859',[0 10],Fi0);
    %Vector de flujos iniciales(con reciclo)
    Fi2=[Fi2(end,1),Fi2(end,2),Fi2(end,3)];       
    
    Flujos_iter=[Fi2];  
    k=1;
    while tol>0.01

        %Ecuacion de reciclo:
        Fa1=Fa0+(R*Fi0(1));
        Fb1=Fb0+(R*Fi0(2));
        Fc1=Fc0+(R*Fi0(3));
        
        %Nuevos Flujos de entrada con reciclo
        Fa0=Fa1;
        Fb0=Fb1;
        Fc0=Fc1;
        
        %Resolucion del TAC por segunda vez
        [Volumen,Fi2]=ode45('SED_1311073_1110859',[0 10],[Fa0,Fb0,Fc0]);
        Fk=[Fi2(end,1),Fi2(end,2),Fi2(end,3)];     %Vector de flujos iter.k(con reciclo)
        Flujos_iter=[Flujos_iter; Fk];             %Tabla de flujos para cada iteracion.    
 
        %Calculo tolerancia
        tol=max(abs(((Flujos_iter(k,:)-Flujos_iter(k+1,:))/Flujos_iter(k+1,:))));
        k=k+1;
        if k>=50
            break
        end
    end
    Convergencia_iteracion=k
    Flujos_iter
    Flujos_salida=Flujos_iter(end,:)

    hold on
    plot([1:k],Flujos_iter(:,1),'*b')
    plot([1:k],Flujos_iter(:,2),'og')
    plot([1:k],Flujos_iter(:,3),'-y')    
end
