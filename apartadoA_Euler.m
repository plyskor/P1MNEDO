S0 = 500;
I0 = 0;
Z0 = 2;
R0 = 0;
a = 0.005;
b = 0.0095;
c = 0.0001;
d =3;

ini=0;
fin=10;
pasos=500;

figure(1);
hold all;

%Metodo de Euler
[S,I,Z,R,t] = euler(@fun,pasos,ini,fin,S0,I0,Z0,R0,a,b,c,d)

sum=S+I+Z+R;
plot(S);
plot(I);
plot(Z);
plot(R);
plot(sum);
title('Solucion numerica');
xlabel('x');
ylabel('y');
legend('Humanos','Infectados','Zombies','Muertos');
hold all;
