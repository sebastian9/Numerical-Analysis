function [X,Y] = RungeKutta
syms x y
fxy = input('Inserte la función f(x,y) igual a dy/dx\n');
x0 = input('Introduzca el valor inicial de x\n');
y0 = input('Introduzca el valor correspondiente de y\n');
xf = input('Introduzca el valor final de x\n');
instruction = ['Elija el enfoque a utilizar\n' ...
    'a) Dividir el intervalo en n pasos de igual tamaño\n' ...
    'b) Utilizar un tamaño de paso\n'];
option = input( instruction, 's');
if option == 'a'
    n = input('Introduzca el número de pasos\n');
    h = (xf-x0)/n;
elseif option == 'b'    
    h = input('Introduzca el tamaño de paso\n');
    n = (xf-x0)/h;
else    
    return
end

i = 1;
X = zeros(1,n);
Y = zeros(1,n);
X(i) = x0;
Y(i) = y0;

for xn = x0:h:xf-h
    yn = Y(i);
    k1 = eval( subs ( subs (fxy,x,xn),y,yn ) );
    k2 = eval( subs ( subs (fxy,x,xn + 0.5*h ),y,yn + 0.5*h*k1 ) );
    k3 = eval( subs ( subs (fxy,x,xn + 0.5*h ),y,yn + 0.5*h*k2 ) );
    k4 = eval( subs ( subs (fxy,x,xn + h ),y,yn + h*k3 ) );
    Y(i+1) = yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    X(i+1) = xn + h;
    i = i + 1;
end

X = transpose(X);
Y = transpose(Y);

fitpoly = fit(X,Y,'poly5')
figure, plot(fitpoly,X,Y);