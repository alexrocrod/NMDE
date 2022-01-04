function [x, iter, resvec, flag] = mygmres(A, b, tol, maxit, x0)
    
    r = b - A * x0;
    v = zeros(maxit + 1, length(r));
    h = zeros(maxit + 1, maxit + 1);
    rho(1) = norm(r,2);
    beta = rho;
    v(1,:) = r / beta; 
    flag = 0;
    b_norm = norm(b);
    thres = tol*b_norm;
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
        rho(k+1) = abs(beta*Q(1,k+1));

        if h(k+1,k) == 0
            flag = -1;
            break
        end
        v(k+1,:) = v(k+1,:) /h(k+1,k);
       
        
    end
    e1 = zeros(k,1);
    e1(1) = 1;
    be1 = beta*e1;
    if flag == -1
        y = H\be1;
    else
        y = R(1:k,1:k)\(Q(1:k,1:k)'*be1);
    end
    x = x0 + v(1:k,:)'*y;

    resvec = rho; 
    iter = k;
        
end
% 
% 
%         k = k + 1
%         v(:,k+1) = A * v(:,k);
%         for j = 1:k
%             h(j,k) = v(:,k+1)' * v(:,j);
%             v(:,k+1) = v(:,k+1) - h(j,k) * v(:,j);
%         end
% 
%         h(k+1,k) = norm(v(:,k+1));
%         
% %         H = h(1:k+1, 1:k);
%         v(:,k+1) = v(:,k+1) / h(k+1,k);         
%        
%     end
%     e0 = zeros(k,1);
%     e1 = e0;
%     e1(1) = beta;
%     if flag == -1
%         y = h(1:k-1,:)\e1;
%     else
%         y = h(1:k,1:k) \ rho(1:k);
%     end
%     x = x0 + v(:,1:k) * y;
% end
