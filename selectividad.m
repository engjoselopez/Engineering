%Calculo de la selectividad de D frente a E

function S=selectividad(k1,k2,k3,fpi,vo)

FA=fpi(1,1);FB=fpi(2,1);FC=fpi(3,1);FD=fpi(4,1);FE=fpi(5,1);FF=fpi(6,1);
CA=FA/vo;CB=FB/vo;CC=FC/vo;CD=FD/vo;CE=FE/vo;CF=FF/vo;
 
rd1=k1*CA*(CB^2);
re2=k2*CA*CD;
rf3=k3*CB*(CC^2);
 
rd=rd1-2*re2+rf3;
re=re2;
  
S=rd/re;  
end function