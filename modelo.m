%Sistema de EDO's
%Aquí está el sistema de ecuaciones diferenciales del sistema dinámico 
%que modela nuestro universo. Las cuatro ecuaciones representan
%respectivamente la variación de las cuatro poblaciones (Humanos Sanos,
%Humanos Infectados, Zombies, Muertos) con respecto de la variable tiempo.

function [dS,dI,dZ,dR]=modelo(S,I,Z,R,a,b,c,d)
    dS=-b*S*Z;          %Poblacion sana
    dI=b*S*Z-d*I;       %Población infectada
    dZ=d*I+c*R-a*S*Z;   %Población zombie
    dR=a*S*Z-c*R;       %Población muerta susceptible de zombificación.
end

