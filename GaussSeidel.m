function [X] = GaussSeidel(A,B)

% new
X = zeros(size(B));
% old
U = zeros(size(B));

flag = 1;
while flag
    for i = 1:length(B)
        sum = 0;
        for j = 1:i-1
            sum = sum - A(i,j) * X(j,1);
        end
        for j = i+1:length(B)
            sum = sum - A(i,j) * U(j,1);
        end
        X(i,1) = ( sum + B(i,1) ) / A(i,i);
    end
    if max( abs( (X-U) / X ) ) < 0.00001
        flag = 0;
    end
    U = X;
end