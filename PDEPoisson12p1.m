% Poisson Equation Finite-Difference; Burden A12.1
% To approximate the solution to the Poisson equation
% ∂2u
% ∂x2 (x, y) +
% ∂2u
% ∂y2 (x, y) = f (x, y), a ≤ x ≤ b, c ≤ y ≤ d,
% subject to the boundary conditions
% u(x, y) = g(x, y) if x = a or x = b and c ≤ y ≤ d
% and
% u(x, y) = g(x, y) if y = c or y = d and a ≤ x ≤ b :

% INPUT endpoints a, b, c, d; integers m ≥ 3, n ≥ 3; tolerance TOL; maximum number of
% iterations N.
% OUTPUT approximations wi,j to u(xi, yj) for each i=1, ..., n−1 and for each j = 1, ... ,
% m − 1 or a message that the maximum number of iterations was exceeded.
function [R] = PDEPoisson12p1(fx,gx,a,b,c,d,m,n,N)

function [r] = fof(a1,a2)
    r = eval( subs ( subs ( fx, x, a1) , y, a2 ) );
    return
end

function [r] = gof(a1,a2)
    r = eval( subs ( subs ( gx, x, a1) , y, a2 ) );
    return
end


% Set h = (b − a)/n;
% k = (d − c)/m.
h = (b-a)/n;
k = (d-c)/m;

% For i = 1, ... , n − 1 set xi = a + ih. (Steps 2 and 3 construct mesh points.)
X = zeros(n-1);
for i = 1:n-1
    X(i,1) = a + i*h;
end

% For j = 1, ... , m − 1 set yj = c + jk
Y = zeros(m-1);
for j = 1:m-1
    Y(j,1)
end

% For i = 1, ... , n − 1
%   for j = 1, ... , m − 1 set wi,j = 0.
W = zeros(n-1,m-1);

% Set λ = h*h/k*k;
% μ = 2(1 + λ);
% l = 1.
lambda = (h*h)/(k*k);
micro = 2*(1 + lambda);
l = 1;

% While l ≤ N do Steps 7–20. (Steps 7–20 perform Gauss-Seidel iterations.)
while l < N
    % Set z = EXPRESSION
    % NORM = |z − w1,m−1|;
    % w1,m−1 = z.
    z = ( - h*h * fof(X(1,1),Y(end)) + gof(a,Y(end) + lambda * gof(X(1,1),d) + lambda * W(1,m-2) + W(2,m-1) )/micro;
    NORM = (z - W(1,m-1));
    W(1,m-1) = z;
    
%     For i = 2, ... , n − 2
%         set z = 
%         − h2f (xi, ym−1) + λg(xi, d) + wi−1,m−1
%         +wi+1,m−1 + λwi,m−2
% 
%         /μ;
%         if |wi,m−1 − z| > NORM then set NORM = |wi,m−1 − z|;
%         set wi,m−1 = z.
    for i = 2:n-2
        z = ( - h*h * fof(X(i,1),Y(m-1,1)) + lambda * gof(X(i,1), d) + W(i-1,m-1) + W(i+1,m-1) + lambda * W(i,m-2) ) / micro;
        if norm( W(i,m-1) - z ) > NORM
            NORM = norm( W(i,m-1) - z );
        end
        W(i,m-1) = z;
        z = ( - h*h * fof(X(n-1,1),Y(m-1,1)) + gof(b,Y(m-1,1)) + lambda * gof(X(n-1,1),d) + W(n-2,m-1) + lambda * W(n-1,m-2) ) / micro;
        if norm( W(n-1,m-1) - z ) > NORM
            NORM = norm( W(n-1,m-1) - z );
        end
        W(n-1,m-1) = z;
        for j = m-2:2
            z = ( - h*h * fof(X(1,1),Y(j,1)) + gof(a,Y(j,1)) + lambda * W(1,j+1) + lambda * W(1,j-1) + W(2,j) ) / micro;
            if norm( W(1,j) - z ) > NORM
                NORM = norm( W(1,j) - z );
            end
            W(1,j) = z;
            for i2 = 2:n-2
                z = ( - h*h * fof(X(i2,1),Y(j,1)) + W(i2-1,j) + lambda * W(i2,j+1) + W(i2+1,j) + lambda * W(i2,j-1) ) / micro;
                if norm( W(i2,j) - z ) > NORM
                    NORM = norm( W(i2,j) - z );
                end
                W(i2,j) = z;
            end
            z = ( - h*h * fof(X(n-1,1),Y(j,1)) + gof(b,Y(j,1)) + W(n-2,j) + lambda * W(n-1,j+1) + lamda * W(n-1,j-1) ) / micro;
            if norm( W(n-1,j) - z ) > NORM
                NORM = norm( W(n-1,j) - z );
            end
            W(n-1,j) = z;
        end
        z = ( - h*h * fof(X(1,1),Y(1,1)) + gof(a,1) + lambda * gof(1,c) + lambda * W(1,2) + W(2,1) ) / micro;
        if norm( W(1,1) - z ) > NORM
            NORM = norm( W(1,1) - z );
        end
        W(1,1) = z;
        for i2 = 2:n-2
            z = ( - h*h * fof(X(i2,1),Y(1,1)) + lambda * gof(X(i2,1),c) + W(i2-1,1) + lambda * W(i2,2) + W(i2+1,1) ) / micro;
            if norm( W(i2,1) - z ) > NORM
                NORM = norm( W(i2,1) - z );
            end
            W(i2,1) = z;
        end
        z = ( - h*h * fof(X(n-1,1),Y(1,1)) + gof(b,Y(1,1)) + lambda * gof(X(n-1,1),c) + W(n-2,1) + lamda * W(n-1,2) ) / micro;
        if norm( W(n-1,1) - z ) > NORM
            NORM = norm( W(n-1,1) - z );
        end
        W(n-1,1) = z;
        if NORM <= N
           R = W;
        end
        l = l+1;
    end
end

end
