%Funcion con las ED para el FPI

function F=FPI(Flujo,V,vo)

k1=0.25; k2=0.10;k3=5;
FA=Flujo(1,1);FB=Flujo(2,1);FC=Flujo(3,1);FD=Flujo(4,1);FE=Flujo(5,1);FF=Flujo(6,1);
CA=FA/vo;CB=FB/vo;CC=FC/vo;CD=FD/vo;CE=FE/vo;CF=FF/vo;

rd1=k1*CA*(CB^2);
re2=k2*CA*CD;
rf3=k3*CB*(CC^2);

ra=-rd1-3*re2;
rb=-2*rd1-rf3;
rc=rd1+re2-2*rf3;
rd=rd1-2*re2+rf3;
re=re2;
rf=rf3;

F(1,1)=ra;
F(2,1)=rb;
F(3,1)=rc;
F(4,1)=rd;
F(5,1)=re;
F(6,1)=rf; 
    
endfunction