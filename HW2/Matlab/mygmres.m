function [x, iter, resvec, flag] = mygmres(A, b, tol, maxit, x0)
    k = 1;
    r = b - A * x0;
    v = zeros(maxit + 1, length(r));
    rho(k) = norm(r,2);
    beta = rho;
    v(k,:) = r / beta; 
    flag = 0;

    while (rho(k) > tol*norm(b,2) & k < maxit) 
        
        v(k+1,:) = A * v(k,:)';
        for j = 1:k
            h(j,k) = v(k+1,:) * v(j,:)';
            v(k+1,:) = v(k+1,:) - h(j,k) * v(j,:);
        end

        h(k+1,k) = norm(v(k+1),2);
        H = h(1:k+1, 1:k);
        [Q,R] = qr(H);
        rho(k+1) = norm(beta*Q(1,k+1));

        if h(k+1,k) == 0
            flag = -1;
            break
        end
        v(k+1,:) = v(k+1,:) / h(k+1,k);
       
        k = k + 1;
    end
    e0 = zeros(k,1);
    e1 = e0;
    e1(1) = beta;
    if flag == -1
        y = H\e1;
    else
        y = R\(beta*Q'*e1)
    end
    x = x0 + v(1:k-1,:)'*y;

    resvec = rho; 
    iter = k-1;
        
end
