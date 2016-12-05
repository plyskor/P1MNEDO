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
d = 3;
%Queremos que haya 1000 nodos y 10 d�as
inicio=0;
dias=10;
pasos=1000;

%Llamamos al m�todo de Adams-Moulton pasando los argumentos en orden, pero al final
%de cada dia...
for i=0:dias-1
   %Vamos avanzando de nodo cada iteracion 
   nodo=i*floor(pasos/dias)+1;
   nextnodo=(i+1)*floor(pasos/dias)+1;
   %Aqu� se aplica A-M
   [S(nodo:nextnodo),I(nodo:nextnodo),Z(nodo:nextnodo),R(nodo:nextnodo),t(nodo:nextnodo)] = adamsMoulton(@modelo,floor(pasos/dias),inicio+i,inicio+i+1,S0,I0,Z0,R0,a,b,c,d);
   %... Aqu� sale el sol y hay pelea
   I0=I(nextnodo)+S(nextnodo)/9;%Infectamos a un noveno de los sanos
   S0=S(nextnodo)-S(nextnodo)/9;%y los quitamos del conteo de sanos
   Z0=Z(nextnodo)/3; %Matamos a dos tercios de zombies
   R0=R(nextnodo)+2*Z(nextnodo)/3; %a�adimos los zombies muertos al conteo de muertos
end
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
title('Apartado 2.b: M�todo de Adams-Moulton (2 pasos)');
xlabel('Nodos');
ylabel('Cantidad');
legend('Humanos sanos','Zombies','Suma');
hold all;
