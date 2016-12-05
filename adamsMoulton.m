%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCION ADAMS-MOULTON 2 PASOS
%
% u = f(t,u)
% u(t0) = u0
% pasos : numero de pasos en los que dividimos el intervalo de integracion
% inicio : tiempo inicial
% T : tiempo final
% S0: valor inicial de S
% Z0: valor inicial de Z
% R0: valor inicial de muertos
% a,b,c: ctes
%
% Variables de Salida:
% u: matriz de lenght(u0) x N que contiene la solucion
% t: vector de tiempos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,I,Z,R,tam] = adamsMoulton(f, pasos , inicio , fin , S0, I0, Z0, R0, a, b, c, d)
	h = (fin-inicio)/pasos;
	%defino tam (el valor del paso) e inicializo las variables donde voy a guardarlos
	tam = [inicio:h:fin];
	S = zeros(1,length(tam));
    I = zeros(1,length(tam));
	Z = zeros(1,length(tam));
	R = zeros(1,length(tam));
	dszr = zeros(4,4); %funcion auxiliar para guardar paso n, n+1 y n+2 de S, Z y R
	szr2aux = zeros(4,1); 
	szr2 = zeros(4,1);

	[S(1:2),I(1:2),Z(1:2),R(1:2),tam(1:2)]=euler(f,1,inicio,inicio+h,S0,I0,Z0,R0,a,b,c,d) %calculamos paso 0 y 1 y los guardamos en las finales

	max_iter= 1000 %stop condition
	tolerancia=0.001; %si la differencia es esta nos contentamos y paramos

	for paso = 1:pasos-1
		[dszr(1,1),dszr(1,2),dszr(1,3),dszr(1,4)]=f(S(paso),I(paso),Z(paso),R(paso),a,b,c,d); %paso n
		[dszr(2,1),dszr(2,2),dszr(2,3),dszr(2,4)]=f(S(paso+1),I(paso+1),Z(paso+1),R(paso+1),a,b,c,d); %paso n+1

		szr2 =[S(paso+1),I(paso+1),Z(paso+1),R(paso+1)] %valores paso n+1

		for i = 1:max_iter %este bucle se hace hasta que la diferencia es pequeña o hasta que nos cansamos (max_iter)
			
			szr2aux=szr2; %guardamos para comparar

			[dszr(3,1),dszr(3,2),dszr(3,3),dszr(3,4)]= f(szr2(1),szr2(2),szr2(3),szr2(4),a,b,c,d); %paso n+2

			szr2=[S(paso+1),I(paso+1),Z(paso+1),R(paso+1)]+h/12*(5*dszr(3,:)+8*dszr(2,:)-dszr(1,:)); %aplicamos al joven Adams

			if (abs(norm(szr2-szr2aux,1))<tolerancia) %%Si la diferencia entre los dos ultimos pasos es pequeña, ya paramos
				break;
			end
		end
		S(paso+2)=szr2(1);
        I(paso+2)=szr2(2);
		Z(paso+2)=szr2(3);
		R(paso+2)=szr2(4);   
	end
end
