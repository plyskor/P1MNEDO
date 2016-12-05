function [u,t]=euler(f,N,t0,T,u0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta funci�n resuelve el problema de valor inicial
% u�=f(t,u)
% u(t0)=u0
% utilizando el m�todo de Euler
%
% [u,t]=euler(N,t0,T,u0)
%
% Variables de Entrada:
%
% f: vector columna. funci�n que rige el sistema de EDO,
% tiene dos argumentos f(t,u) donde t es escalar
% y u vector columna.
% N: n�mero de pasos en los que dividimos el intervalo de
% integraci�n
% t0: tiempo inicial
% T: tiempo final
% u0: vector columna. Dato inicial
%
% Variables de Salida:
% u: matriz de length(u0) x N que contiene la soluci�n
% t: vector de tiempos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defino el paso
h=(T-t0)/N;
% defino u y t
t=[t0:h:T];
u=zeros(length(u0),length(t));
% identifico el dato inicial
u(:,1)=u0;
% aplicamos el algoritmo de Euler
for n=1:N
    u(:,n+1)=u(:,n)+h*f(t(n),u(:,n));
end
end


% Defino la funci�n f de mi PVI yprima=f(x,y(x))

function y=f(t,v)
% v: vector columna de tama�o 2
% t: escalar
% y: vector columna de tama�o 2
y=[-2,1;1,-2]*v+[2*sin(t); 2*(cos(t)-sin(t))];
end

function U= Iteracion_Funcional(f,tn,un,h,A,c,tol,itmax)
    s=length(c);
    e=ones(s,1);
    % El siguiente comando sirve a para dar el candidato inicial. Crea una
    % matriz con s columnas y 2 filas. Dados dos vectores A y B kron(A,B)
    % es el producto tensorial de A y B. Es decir, la matriz C, dada por el
    % producto tensorial de A y B tiene como entradas C_{ij}= A_i B_j. Hay
    % diferentes maneras de crear esta matriz. Por ejemplo podriamos
    % utililzar un comando for 
    % for i=1:s
    %for j=1:2
    %U(j,i)= un(j).
    %end
    %end
    %tambi�n podr�amos haber eligido U=zeros(length(un),s), que es una
    %matriz de ceros.
    U=kron(un,e);
    %iniciamos un contador de iteraciones
    it=0;
    %iniciamos el parametro norma.
    norma=tol+1; %Nos garantizamos hacer por lo menos una iteraci�n
    while (norma>tol && it< itmax) 
        it=it+1;
        auxU=U;
        for i=1:s
        for j=1:s 
        V=A(i,j)*U(:,j);
        end
        U(:,i)=f(tn+c(i)*h,un+h*V);
        end
        norma=norm(U-auxU)/norm(U);
    end
end

% Defino los parametros  iniciales.


% Defino el dato inicial. Ser� un vector columna.
u0=[2; 3];
%Tiempo inicial
t0=0;
%Tiempo final
tf=10;
% N�mero de nodos.
N=500;

%Defino el paso h
h=(tf-t0)/N;



% Llamo a la funci�n que aplica a mi sistema el m�todo de Euler. Esta
% funci�n devuelve la solucion u y el vector de tiempos t.
 [u,t]=euler(@f,N,t0,tf,u0);
 
 %Dibujo la solucion num�rica obtenida.
 figure(1);hold all;
 plot(u(1,:), u(2,:)); hold all

 

 
%Dibujo la soluci�n exacta

uexac=solexac_aPVI(t);
figure(1) 
hold all
plot(uexac(1,:), uexac(2,:),'r');

%Representamos mediante gr�ficas doblemente logaritmicas las poligonales
%(1/h,error);

for N=500:50:1000
    i=i+1;
    h((N-500)/50+1)=(tf-t0)/N;
    %Euler
    [u,t]=euler(@f,N,t0,tf,u0);
    uexac=solexac_aPVI(t); %solucion exacta
    error(1,(N-500)/50+1)=norm((u-uexac),1); 
end

% calculo la pendiente de la recta logaritmo del error= m*logaritmo de h +
% constante
pendiente(1)=mean(diff(log(error(1,:)))./diff(log(h)));

%Dibujo el logaritmo del error frente al logaritmo de h
figure(2); loglog(h,error); title (['La pendiente es =',num2str(pendiente)]);

% Defino los parametros iniciales

u0=[2; 3];
t0=0;
tf=10;
N=500;h=(tf-t0)/N;

% Defino los parametros del Runge-kutta

c=[0,1/2,1/2,1];
b=[1/6,1/3,1/3,1/6];
A=[0,0,0,0;1/2,0,0,0;0,1/2,0,0;0,0,1,0];

hold on
%M�todo de runge-kutta
[u,t]=RKexplicito(@f,tf,t0,N,u0,b,c,A);
figure(1);hold all;
plot(u(1,:), u(2,:)); hold all

%soluci�n exacta

uexac=solexac_aPVI(t);
figure(1);hold all;
plot(uexac(1,:), uexac(2,:));

%Representamos mediante gr�ficas doblemente logaritmicas las poligonales
%(1/h,error);

for N=500:50:1000
    i=i+1;
    h((N-500)/50+1)=(tf-t0)/N;
    %Euler
    [u,t]=RKexplicito(@f,tf,t0,N,u0,b,c,A);
    uexac=solexac_aPVI(t); %solucion exacta
    error(1,(N-500)/50+1)=norm((u-uexac),1); 
end
pendiente(1)=mean(diff(log(error(1,:)))./diff(log(h)));

figure(2); loglog(h,error); title (['La pendiente es =',num2str(pendiente)]);



% Definolos parametros iniciales
u0=[2; 3];
t0=0;
tf=10;
h=0.01; 
itmax=100;
tol=1d-10;


%Aqu� os dejo los comandos para ejecutar 4 runge-kuttas distintos

%Metodo de Gauss: Regla impl�cita del punto medio
A=1/2; b=1; c=1/2;

[u,t]=RKimplicito(@f,tf,t0,h,u0,A,b,c,tol,itmax);
figure(1); plot(u(1,:),u(2,:),'-'); hold all


%Euler Impl�cito
% A=1; b=1; c=1;
% 
% [u,t]=RKimplicito(@f,tf,t0,h,u0,A,b,c,tol,itmax);
% figure(1); plot(u(1,:),u(2,:),'-'); hold all


%M�todo de Gauss
% c=[1/2-sqrt(15)/10; 1/2; 1/2+sqrt(15)/10];
% b=[5/18, 4/9, 5/18];
% A=[5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30; 5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24; 5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
% 
% [u,t]=RKimplicito(@f,tf,t0,h,u0,A,b,c,tol,itmax);
% figure(1); plot(u(1,:),u(2,:)); hold all


%Lobatto IIIA. Regla Trapezoidal
% A=[0,0; 1/2,1/2];
% b=[1/2,1/2];
% c=[0;1];
% 
% [u,t]=RKimplicito(@f,tf,t0,h,u0,A,b,c,tol,itmax);
% figure(1); plot(u(1,:),u(2,:),'--'); hold all



% %Representar la solucion aproximada junto con la soluci�n exacta.

 uexac=solexac_aPVI(t);
 figure(1); plot(uexac(1,:), uexac(2,:)); hold all
 

 
 function [u,t]=RKexplicito(f,tf,t0,N,u0,b,c,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta funci�n resuelve el problema de valor inicial
%       u=f(t,u)
%       u(t0)=u0
% utilizandodo un m�todo Runge-Kutta explicito
%
%           [u,t]=RKexplicito(f,tf,t0,N,u0,b,c,A)
%
%  Variables de Entrada:
%
%       f: vector columna. funci�n que rige el sistema de EDO, 
%          tiene dos argumentos f(t,u) donde t es escalar
%          y u vector columna.   
%       N: n�mero de pasos en los que dividimos el intervalo de
%          integraci�n
%       t0: tiempo inicial 
%       tf: tiempo final
%       u0: vector columna. Dato inicial
%       b,c,A: tablero de BUTCHER. 
%           A: matriz cuadrada s
%           c: vector de tama�o s
%           b: vector de tama�o s
%
%  Variables de Salida:
%
%       u: matriz de length(u0) x length(t) que contiene la soluci�n
%       t: vector de tiempos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=(tf-t0)/N;
t=[t0:h:tf];
u=zeros(length(u0),length(t));
u(:,1)=u0; %En la columna 1 guardamos el dato  inicial 
%
%Aplicamos el algoritmo Runge-Kutta
s=length(b); %n�mero de ks
k=zeros(length(u0),s);

for n=1:N
    k(:,1)=f(t(n)+h*c(1),u(:,n));
    for i=2:s
        k(:,i)=f(t(n)+c(i)*h,u(:,n)+h*k(:,1:i-1)*transpose(A(i,1:i-1)));
    end
    u(:,n+1)=u(:,n)+h*k*b;
end

end
 
 function [u,t]=RKimplicito(f,tf,t0,h,u0,A,b,c,tol,itmax)
%[u,t,it]=RKimplicito(f,df,tf,t0,h,u0,A,b,c,Met,tol,itmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Esta funci�n resuelve el problema de valor inicial
%       u=f(t,u)
%       u(t0)=u0
% utilizandodo un m�todo Runge-Kutta implicito
%
%          
%
%  Variables de Entrada:
%
%       f: vector columna. funci�n que rige el sistema de EDO, 
%          tiene dos argumentos f(t,u) donde t es escalar
%          y u vector columna.   
%       h : tama�o del paso 
%       t0: tiempo inicial 
%       tf: tiempo final
%       u0: vector columna. Dato inicial
%       b,c,A: coeficientes del tablero de BUTCHER.       
%       
%       tol: tolerancia para las iteraciones
%       itmax: Numero m�ximo de iteraciones
%       
%  Variables de Salida:
%
%       u: matriz de length(u0) x length(t) que contiene la soluci�n
%       t: vector de tiempos
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=[t0:h:tf]; N=length(t);
m=length(u0);
u=zeros(m,N);
u(:,1)=u0; %En la columna 1 guardamos la soluci�n inicial 


s=size(A,1);


for n=1:N-1        
U=Iteracion_Funcional(f,t(n),u(:,n),h,A,c,tol,itmax);
u(:,n+1)=u(:,n)+h*U*transpose(b); 
end
end

% Defino la soluci�n exacta que para el caso particular de la funci�n f
% definida en el archivo f.m conozco.
function uexac=solexac_aPVI(t)
    uexac(1,:)=  2*exp(-t)+sin(t);
    uexac(2,:) =2*exp(-t)+cos(t) ;
end
