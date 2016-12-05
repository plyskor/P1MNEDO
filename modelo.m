%Sistema de EDO's
%Aqu� est� el sistema de ecuaciones diferenciales del sistema din�mico 
%que modela nuestro universo. Las cuatro ecuaciones representan
%respectivamente la variaci�n de las cuatro poblaciones (Humanos Sanos,
%Humanos Infectados, Zombies, Muertos) con respecto de la variable tiempo.

function [dS,dI,dZ,dR]=modelo(S,I,Z,R,a,b,c,d)
    dS=-b*S*Z;          %Poblacion sana
    dI=b*S*Z-d*I;       %Poblaci�n infectada
    dZ=d*I+c*R-a*S*Z;   %Poblaci�n zombie
    dR=a*S*Z-c*R;       %Poblaci�n muerta susceptible de zombificaci�n.
end

