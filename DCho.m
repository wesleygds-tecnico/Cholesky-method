function [] = DCho(A)
n = length(A);
B = eig(A);
%{
for k = 1:n
    a = A(1:k,1:k);
    if ((B(k,1)>0)&&(det(a)>0))
        a = 1;
    else 
        fprintf('A matriz não é SP \n')
        a = 0;
        break
    end
end
%}
a = 1;
if (a == 1)
    fprintf('A matriz é SP \n')
    L = zeros(n);
    for i = 1:n 
        L(i,i) = sqrt(A(i,i)-(L(i,:))*L(i,:)');
        for k = (i+1):n
            L(k,i) = (A(k,i)-L(i,:)*L(k,:)')/L(i,i);
        end
    end
    disp(L)
    disp(L')
    C = L*L';
    fprintf('Verifiação \n')
    disp(C)
end      
%disp(G')
%C = G'*G;
%disp(C)
end