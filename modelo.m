%Sistema de EDO's
%En este documento queda definido el sistema de Ecuaciones Diferenciales
%que nos sirve como modelo matemático (sistema dinámico) de las poblaciones
%que queremos estudiar.
function [dS,dI,dZ,dR]=fun(S,I,Z,R,a,b,c,d)
    dS=-b*S*Z;          %Poblacion Sana
    dI=b*S*Z-d*I;       %Población Infectada
    dZ=d*I+c*R-a*S*Z;   %Población Zombie
    dR=a*S*Z-c*R;       %Población Muerta Zombificable
end

