%Introducimos los valores iniciales del PVI que aparecen en el enunciado
%para las poblaciones
S0 = 500;
I0 = 0;
Z0 = 2;
R0 = 0;
%y para los parametros de crecimiento/extincion
a = 0.005;
b = 0.0095;
c = 0.0001;
d =3;
%Queremos que haya 1000 nodos y 10 días
ini=0;
fin=10;
pasos=1000;

%Llamamos al método de Euler pasando los argumentos en orden
[S,I,Z,R,t] = euler(@modelo,pasos,ini,fin,S0,I0,Z0,R0,a,b,c,d)

%Construimos una gráfica, donde mostraremos las dos poblaciones que
%nos piden (S'z y Z's) y además también representaremos el valor sum,
%que es el cómputo total de individuos
figure(1);
hold all;
sum=S+I+Z+R;
plot(S);
plot(Z);
plot(sum);

%Etiquetamos la gráfica y listo
title('Apartado 1.a: Método de Euler');
xlabel('Nodos');
ylabel('Cantidad');
legend('Humanos sanos','Zombies','Suma');
hold all;
