function [x, iter, resvec, flag] = mygmres(A, b, tol, maxit, x0)
    r = b - A * x0;
    v = zeros(maxit + 1, length(r));
    h = zeros(maxit + 1, maxit + 1);
    rho = zeros(maxit+1,1);

    beta = norm(r,2);
    rho(1) = beta;
    v(1,:) = r / beta;
    thres = tol * norm(b);
    
    flag = 0;
    k = 0;
    while (rho(k+1) > thres && k < maxit) 
        k = k + 1;
        v(k+1,:) = A * v(k,:)';
        for j = 1:k
            h(j,k) = v(k+1,:) * v(j,:)';
            v(k+1,:) = v(k+1,:) - h(j,k) * v(j,:);
        end

        h(k+1,k) = norm(v(k+1,:));
        H = h(1:k+1, 1:k);
        [Q,R] = qr(H);
        rho(k+1) = abs(beta * Q(1,k+1));

        if h(k+1,k) == 0
            flag = -1;
            break
        end
        v(k+1,:) = v(k+1,:) / h(k+1,k);
       
        
    end
    e1 = zeros(k,1); e1(1) = 1; be1 = beta * e1;
    if flag == -1
        y = H\be1;
    else
        y = R(1:k,1:k) \ (Q(1:k,1:k)' * be1);
    end
    x = x0 + v(1:k,:)' * y;

    resvec = rho(1:k+1); 
    iter = k;
        
end
