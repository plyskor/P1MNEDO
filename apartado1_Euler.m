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
%Queremos que haya 1000 nodos y 10 d�as
ini=0;
fin=10;
pasos=1000;

%Llamamos al m�todo de Euler pasando los argumentos en orden
[S,I,Z,R,t] = euler(@modelo,pasos,ini,fin,S0,I0,Z0,R0,a,b,c,d)

%Construimos una gr�fica, donde mostraremos las dos poblaciones que
%nos piden (S'z y Z's) y adem�s tambi�n representaremos el valor sum,
%que es el c�mputo total de individuos
figure(1);
hold all;
sum=S+I+Z+R;
plot(S);
plot(Z);
plot(sum);

%Etiquetamos la gr�fica y listo
title('Apartado 1.a: M�todo de Euler');
xlabel('Nodos');
ylabel('Cantidad');
legend('Humanos sanos','Zombies','Suma');
hold all;
