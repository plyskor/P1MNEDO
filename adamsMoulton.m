function [S,I,Z,R,VectorTiempoDiscreto] = adamsMoulton(funcion, pasos , start , stop , S0, I0, Z0, R0, a, b, c, d)
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
    %Construimos la matriz auxiliar para guardar los pasos n, n+1 y n+2 de
    %las cuatro poblaciones
	matrizPasos = zeros(3,4); 
    %Los vectores fila para la iteracion en cada paso y para la comparación
    %para comprobar que hemos afinado bastante (tolerancia)
	prevEnemasuno = zeros(4,1); 
	enemasuno = zeros(4,1);
    %Usamos Euler para los pasos 0 y 1 y lo guardamos.
	[S(1:2),I(1:2),Z(1:2),R(1:2),VectorTiempoDiscreto(1:2)]=euler(funcion,1,start,start+h,S0,I0,Z0,R0,a,b,c,d);
    %Definimos las condiciones de parada
	max_iter = 10000 %Para evitar un bucle infinito
	tol=1.0e+005; %y para dejar de contar si ya estamos muy bien aproximados
    %Comenzamos el bucle para el método de Adams-Moulton de dos pasos
	for paso = 1:pasos-1
        %Hacemos el paso n y el n+1
		[matrizPasos(1,1),matrizPasos(1,2),matrizPasos(1,3),matrizPasos(1,4)]=funcion(S(paso),I(paso),Z(paso),R(paso),a,b,c,d); 
		[matrizPasos(2,1),matrizPasos(2,2),matrizPasos(2,3),matrizPasos(2,4)]=funcion(S(paso+1),I(paso+1),Z(paso+1),R(paso+1),a,b,c,d); 
        %Metemos los valores del paso n+1 en su vector fila
		enemasuno =[S(paso+1),I(paso+1),Z(paso+1),R(paso+1)] 
        %Iteramos Adams-Moulton hasta que se cumpla una condición de parada
		for i = 1:max_iter 
			prevEnemasuno=enemasuno; 
            %hacemos el paso n+2
			[matrizPasos(3,1),matrizPasos(3,2),matrizPasos(3,3),matrizPasos(3,4)]= funcion(enemasuno(1),enemasuno(2),enemasuno(3),enemasuno(4),a,b,c,d);
			enemasuno=[S(paso+1),I(paso+1),Z(paso+1),R(paso+1)]+h/12*(5*matrizPasos(3,:)+8*matrizPasos(2,:)-matrizPasos(1,:)); 
            %Comprobamos que no hemos sido ya suficientemente finos
			if (abs(norm(enemasuno-prevEnemasuno,1))<tol);
				%si lo hemos sido paramos, si no, vuelta
                break;
			end
        end
        %vamos metiendo los valores en los vectores solucion
		S(paso+2)=enemasuno(1);
        I(paso+2)=enemasuno(2);
		Z(paso+2)=enemasuno(3);
		R(paso+2)=enemasuno(4);   
	end
end
