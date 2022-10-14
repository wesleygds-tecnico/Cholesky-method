function [] = L7Ex3()
C = [1 0 1;1 1 0;1 1 0];
b = [2;2;2];
A = C'*C;
D = C'*b;
n = length(A);
B = eig(A)
for k  = 1:n
    a = A(1:k,1:k);
    if (((B(k,1)>=0) && (det(a))>0) && (((conj(B(k,1)) + B(k,1))/2) == B(k,1)))
        b = 1;
    else
        b = 0;
        fprintf('Não é possível Cholesky decompor a matriz \n')
        return
    end
end
    
%1° encontrar elemento 1,1
%2° encontrar elementos da coluna 1
%3° encontrar elemento 2,2
%4° encontrar elementos da coluna 2....
if (b == 1)
    fprintf('É possível Cholesky decompor a matriz \n')
    G = zeros(n);
    G(1,1) = (A(1,1))^(1/2); 
    G(2:n,1) = A(2:n,1)/G(1,1);
    for k = 2:(n) %coluna
        soma = 0;
        soma1 = 0;
        for i = (k-1):-1:1
            soma1 = soma1 + G(k,i)*G(k,i);
        end
        G(k,k) = (A(k,k)- soma1)^(1/2);
        for i = (n-1):-1:2
            soma = soma + G(i+1,k-1)*G(i,k-1);
        end
        for i = (k+1):n
            G(i,k) = (A(i,k)-soma)/(G(k,k));
        end
    end
    fprintf('Resultado \n')
    disp(G)
    c = G*G';
    fprintf('Verificação \n')
    disp(c)
end

y = zeros(n,1); 
y(1,1) = D(1,1)/G(1,1);
for i = 2:(n)
    soma = 0;
    for j = 1:i-1
        soma = soma + (y(j)*G(i,j));
    end
    y(i) = (D(i)-soma)/G(i,i);
end
U = G';
x = zeros(n,1); 
x(n,1) = y(n,1)/U(n,n);
for i = n-1:-1:1
    soma = 0;
    for j = (i+1):n
        soma = soma + (x(j)*U(i,j));
    end
    x(i) = (y(i)-soma)/(U(i,i));
end
%fprintf('Solução do sistema 1 \n')
%disp(y)
%fprintf('Verificação do sistema 1 \n')
%a = G*y;
%disp(a)
%fprintf('Solução do sistema 2 \n')
%a = G'*x;
%disp(x)
%fprintf('Verificação do sistema 2 \n')
%disp(a)
fprintf('Solução \n')
disp(x)
fprintf('Verificação geral \n')
a = C*x;
disp(a)
end
