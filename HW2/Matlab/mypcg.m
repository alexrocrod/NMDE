function [x, resvec, iter] = mypcg(A, b, tol, maxit, L)
    M1 = L;
    M2 = L';
    M = M1 * M2;
    resvec = zeros(maxit+1, 1);
        
    x = zeros(length(b),1);
    r = b - A * x;
    z = M \ r;
    p = z;
    rsold = r' * z;
    
    for iter = 1:maxit 
        Ap = A * p;
        alpha = rsold / (Ap' * p);
        x = x + alpha * p;
        r = r - alpha * Ap;

        z = M2 \ (M1 \ r);

        rsnew = r' * z;

        resvec(iter) = sqrt(rsnew);
        if sqrt(rsnew) < tol
              break
        end

        p = r + (rsnew / rsold) * p;
        rsold = rsnew; 
    end
    resvec = resvec(1:iter+1);
end
