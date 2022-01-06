function [x, iter, resvec, flag] = myprecgmres(A, b, tol, maxit, x0, L, U)
    LU = L * U;
    r0 = L \ (b - A * x0);
    r = r0;

    h = zeros(maxit + 1, maxit + 1);
    v = zeros(maxit + 1, length(r));
    rho = zeros(maxit + 1, 1);

    beta = norm(r,2);
    rho(1) = beta;
    v(1,:) = r / beta;     
    thres = tol * normest(LU\b);

    flag = 0;
    k = 0;
    while (normest(LU\r) > thres && k < maxit)
        k = k + 1;
        v(k+1,:) = L \ (A * (U \ v(k,:)'));

        for j = 1:k
            h(j,k) = v(k+1,:) * v(j,:)';
            v(k+1,:) = v(k+1,:) - h(j,k) * v(j,:);
        end

        h(k+1,k) = norm(v(k+1,:));
        if h(k+1,k) == 0
            flag = -1;
            break
        end

        H = h(1:k+1, 1:k);
        [Q,R] = qr(H);
        v(k+1,:) = v(k+1,:) / h(k+1,k);

        e1 = zeros(k,1); e1(1) = 1; be1 = beta * e1;

        y = R(1:k,1:k) \ (Q(1:k,1:k)' * be1);
        x = x0 + U \ v(1:k,:)' * y;
        r = L \ (b - A * x);

        rho(k+1) = norm(r);
    end
    if flag == -1
        H = h(1:k+1, 1:k);
        e1 = zeros(k,1); e1(1) = 1; be1 = beta * e1;
        y = H \ be1;
        x = x0 + v(1:k,:)'*y;
    end
    resvec = rho(1:k+1); 
    iter = k;
        
end
