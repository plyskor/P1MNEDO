function [S,I,Z,R,tam] = euler(funcion, pasos, inicio, fin, S0, I0, Z0, R0, a, b, c, d)	
    h = (fin-inicio)/pasos;
	tam = [inicio:h:fin];
	S = zeros(1,length(tam));
    I = zeros(1,length(tam));
	Z = zeros(1,length(tam));
	R = zeros(1,length(tam));
	S(1) = S0;
    I(1) = I0;
	Z(1) = Z0;
	R(1) = R0;
	for paso = 1:pasos
		[prevS, prevI, prevZ, prevR]=funcion(S(paso),I(paso),Z(paso),R(paso),a,b,c,d);
		S(paso+1) = S(paso) + h*prevS;
        I(paso+1) = I(paso) + h*prevI;
		Z(paso+1) = Z(paso) + h*prevZ;
		R(paso+1) = R(paso) + h*prevR;
	end
end
