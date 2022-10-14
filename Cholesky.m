function [] = Cholesky(A)
n = length(A);
B = eig(A);
%{
for k  = 1:n
    a = A(1:k,1:k);
    %if (((B(k,1)>0) && (det(a))>0) && (((conj(B(k,1)) + B(k,1))/2) == B(k,1)))
    if (((det(a))>0) && (((conj(B(k,1)) + B(k,1))/2) == B(k,1)))    
        b = 1;
    else
        b = 0;
        fprintf('Não é possível Cholesky decompor a matriz \n')
        return
    end
end
%}

b = 1
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
    C = G*G';
    fprintf('Verificação \n')
    disp(C)
end
end    