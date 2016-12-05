S0 = 500;
I0 = 0;
Z0 = 2;
R0 = 0;
a = 0.005;
b = 0.0095;
c = 0.0001;
d=3;

t0=0;
dias=10;

pasos=500;

%Metodo de Euler
for i=0:dias-1
   ini=i*floor(pasos/dias)+1;
   fin=(i+1)*floor(pasos/dias)+1;

   [S(ini:fin),I(ini:fin),Z(ini:fin),R(ini:fin),t(ini:fin)] = euler(@fun,floor(pasos/dias),t0+i,t0+i+1,S0,I0,Z0,R0,a,b,c,d);

   S0=S(fin)-S(fin)/9;
   I0=I(fin)+S(fin)/9;
   Z0=Z(fin)/3; %Matamos a dos tercios de zombies
   R0=R(fin)+2*Z(fin)/3; %Los zombis que han muerto se añaden a los muertos.
end

figure(1);
hold all;
n=S+I+Z+R;
plot(S);
plot(I);
plot(Z);
plot(R);
plot(n);
title('Solucion numerica');
xlabel('x');
ylabel('y');
legend('Humanos','Infectados','Zombies','Muertos','Suma');
hold all;
