function [S,I,Z,R,VectorTiempoDiscreto] = euler(funcion, pasos, start, stop, S0, I0, Z0, R0, a, b, c, d)	
    %Calculamos el ancho de paso
    h = (stop-start)/pasos;                 
    %Construimos un vector con el tiempo discretizado en nodos, que sirve tambien
    %para saber el tamaño de los vectores solucion
	VectorTiempoDiscreto = [start:h:stop];
    %Construimos los vectores donde se devuelven los valores
	S = zeros(1,length(VectorTiempoDiscreto));    
    I = zeros(1,length(VectorTiempoDiscreto));      
	Z = zeros(1,length(VectorTiempoDiscreto));      
	R = zeros(1,length(VectorTiempoDiscreto));
    %Introducimos los valores iniciales
	S(1) = S0;        
    I(1) = I0;    
	Z(1) = Z0;     
	R(1) = R0;     
    %y a iterar en cada paso....
	for paso = 1:pasos
		[prevS, prevI, prevZ, prevR]=funcion(S(paso),I(paso),Z(paso),R(paso),a,b,c,d);
		S(paso+1) = S(paso) + h*prevS;
        I(paso+1) = I(paso) + h*prevI;
		Z(paso+1) = Z(paso) + h*prevZ;
		R(paso+1) = R(paso) + h*prevR;
    end
end
