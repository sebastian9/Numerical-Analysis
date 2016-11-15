function [R,fitpoly] = DiferenciasFinitasEDO

syms x y
fprintf('Método por diferencias finitas para el cálculo de PVF de EDO de segundo orden:\n\n');
fprintf('Considere la ecuación de la forma y"" + P(x)y" + Q(x)y = f(x)\n');
fprintf('y condiciones de frontera y(a)=A, y(b) = B\n\n');
px = input('Introduzca el valor de P(x)\n');
qx = input('Introduzca el valor de Q(x)\n');
fx = input('Introduzca el valor de f(x)\n');
a = input('Introduzca el valor de a\n');
ya = input('Introduzca el valor de y(a)\n');
b = input('Introduzca el valor de b\n');
yb = input('Introduzca el valor de y(b)\n');
n = input('A continuación introduzca el número de nodos a incluir en el intervalo\n');
h = (b-a)/n;

%(1 + h*Pi/2)*yip + (-2 + h*h*Qi)*yi + (1 - h*Pi/2)*yim == h*h*fi
A = zeros(n-1);
B = zeros(n-1,1);

for i = 1:n-1
    Ci = -2 + h*h * eval( subs( qx, x, a + h * i ) );
    Cip = 1 + h * eval( subs( px, x, a + h * i ) / 2 );
    Cim = 1 - h * eval( subs( px, x, a + h * i ) / 2 );
    for j = 1:n-1
        if j == i
            A(j,i) = Ci;
        elseif j == i + 1
            A(j,i) = Cim;
        elseif j == i - 1
            A(j,i) = Cip;
        end
    end
    B(i,1) = h*h * eval( subs( fx, x, a + h * i ) );
    if i == 1
        B(i,1) = B(i,1) - Cim * ya;
    elseif i == n-1
        B(n-1,1) = B(i,1) - Cip * yb;
    end
end

Y = [ya; GaussSeidel(A,B); yb];
X = zeros(n,1);
for i = 1:n+1
    X(i,1) = a + h*(i-1);
end

R = [X Y];

fitpoly = fit(X,Y,'poly5');
figure, plot(fitpoly,X,Y);