%Sistema de EDO's
%En este documento queda definido el sistema de Ecuaciones Diferenciales
%que nos sirve como modelo matem�tico (sistema din�mico) de las poblaciones
%que queremos estudiar.
function [dS,dI,dZ,dR]=fun(S,I,Z,R,a,b,c,d)
    dS=-b*S*Z;          %Poblacion Sana
    dI=b*S*Z-d*I;       %Poblaci�n Infectada
    dZ=d*I+c*R-a*S*Z;   %Poblaci�n Zombie
    dR=a*S*Z-c*R;       %Poblaci�n Muerta Zombificable
end

