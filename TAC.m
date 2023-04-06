%Funcion con los balances igualados a cero para el TAC


function F=TAC(C,CAO,CBO,k1,k2,k3,Tau)
CA=C(1);
CB=C(2);
CC=C(3);
CD=C(4);
CE=C(5);
CF=C(6);

rd1=k1*CA*(CB^2);
re2=k2*CA*CD;
rf3=k3*CB*(CC^2);

ra=-rd1-3*re2; 
rb=-2*rd1-rf3;
rc=rd1+re2-2*rf3;
rd=rd1-2*re2+rf3;
re=re2;
rf=rf3;

F(1,1)=CAO-CA+Tau*ra;
F(2,1)=CBO-CB+Tau*rb;
F(3,1)=-CC+Tau*rc;
F(4,1)=-CD+Tau*rd;
F(5,1)=-CE+Tau*re;
F(6,1)=-CF+Tau*rf;

 end