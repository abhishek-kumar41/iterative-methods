function [x_k1, gauss_iter, gauss_conv_check] = Gauss_Siedel_Method(A, b, tolerance, x_k)
m = size(A,1);
x_k1 = zeros(m,1);

for i = 1:m
    sum = 0;
    for j = 1:m
        if j < i
            sum = sum + A(i,j)*x_k1(j,1);
        end
        if j > i
            sum = sum + A(i,j)*x_k(j,1);
        end
    end
    x_k1(i,1) = -sum/A(i,i) + b(i,1)/A(i,i);
end

gauss_iter = 1;
gauss_conv_check = 0;
while (norm(x_k1 - x_k)/norm(x_k1)) >= tolerance 
    x_k = x_k1;
    for i = 1:m
        sum = 0;
        for j = 1:m
            if j < i
                sum = sum + A(i,j)*x_k1(j,1);
            end
            if j > i
                sum = sum + A(i,j)*x_k(j,1);
            end
        end
        x_k1(i,1) = -sum/A(i,i) + b(i,1)/A(i,i);
    end

    gauss_iter = gauss_iter + 1;
    if norm(x_k1) > 1.0e+100
        gauss_conv_check = -1;
        break
    end
end

end

