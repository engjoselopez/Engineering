function F=Balance_TAC(C)

Cao=1.5;
Cbo=2.0;
k1=0.25;
k2=0.10;
k3=5.00;
T=2.5;

r1=k1*C(1)*(C(2)^2);
r2=k2*C(1)*C(4);
r3=k3*C(2)*(C(3)^2);

F(1)=Cao-C(1)-T*(r1+3*r2);
F(2)=Cbo-C(2)-T*(2*r1+r3);
F(3)=-C(3)+T*(r1+r2-2*r3);
F(4)=-C(4)+T*(r1+r3-2*r2);