function [X, Y] = EulerMethod()
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
    m = 0.5 * ( eval( subs (subs (fxy,x,xn),y,yn ) ) + eval( subs (subs (fxy,x,xn),y,yn ) ) );
    Y(i+1) = yn + h * m;
    X(i+1) = xn + h;
    i = i + 1;
end

X = transpose(X);
Y = transpose(Y);

fitpoly = fit(X,Y,'poly9')
figure, plot(fitpoly,X,Y);