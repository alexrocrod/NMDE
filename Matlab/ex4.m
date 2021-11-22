%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all

%% Question 4
y_exact = load("accurate_solution.txt");

method = 0; % 0 - RK4 stability,  1 - ODE45, 2 - Crank-Nicolson, 3 - BDF3

% y' = -Ay;
%% Constants
nx = 100;

G = numgrid ( 'S' , nx ) ;
A = delsq(G) * ( nx-1)^2 ;
lambda = -eigs(A,1,'lm');

T = 0.1;
h=T/length(A);
ts = h:h:T;

N = length(ts);
y = zeros(N,N);
y(1,:) = ones(1,N);

f = @(x) -A*x';

ords = 3:5;
hs1 = 10.^-ords;

%% Stability with RK4
if method == 0
    format shortEng
    h_max_teo = -2.78529/lambda;
    hs = linspace(h_max_teo*0.9, h_max_teo*1.1,10);
%     hs = linspace(h_max_teo*0.99, h_max_teo*1.01,10);
%     hs = linspace(h_max_teo*0.999, h_max_teo*1.001,10);
    errors = zeros(10,1);
    ys = zeros(10,N);
    ih = 1;
    for h=hs
        fprintf('RK4, i=%d,  h=%d\n', ih, h);
        tic
        for i=2:N
            k1 = f(y(i-1,:));
            k2 = f(y(i-1,:)+k1'.*h/2);
            k3 = f(y(i-1,:)+h/2.*k2');
            k4 = f(y(i-1,:)+h.*k3');
                
            y(i,:) = y(i-1,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4)';
        end
        yf = y(end,:)';
        errors(ih) = norm(yf - y_exact, inf);
        ys(ih,:) = y(end,:);
        tab = [hs' errors]
        if any(isnan(ys(ih,:)))
            fprintf('NaN found for h=%d\n',h);
            break
        end
        ih = ih + 1;
        toc

    end
    plot(hs, errors)
    return %end sricpt 
end

%% Solve with ODE45
if method == 1
    disp('Starting ODE45')
    tic
    tspan = [0 0.1];
    [t,y] = ode45(@(t,y) -A*y, tspan, y(1,:));
    yf = y(end,:)';
    error = yf - y_exact;
    fprintf('ODE45, Nsteps = %d, error = %d \n', length(t), norm(error, inf))
    toc
    plot(ts, abs(error))
    return %end sricpt 
end


%% Solve with CN
if method == 2
    disp('Starting CN')

    for h = hs1
        T = 0.1;
        ts = h:h:T;
        N2 = length(ts);
        y = zeros(N2,N);
        y(1,:) = ones(1,N);
        yi = y(1,:)';

        tic
        A2 = eye(N) + A * h/2;
        b = eye(N) - A * h/2;

        fprintf('ConjGrad , h=%d, Nsteps = %d \n', h, N2)
        tol = h^3;
        for i=2:N2
            b2 = b * yi;
            yi = conjgrad(A2, b2, yi, tol);
        end
        error = yi-y_exact;
        fprintf('ConjGrad, Nsteps = %d, error = %d \n', N2, norm(error, inf))
        
        toc   
    end
end


%% Solve with BDF3
% y(i+3) - 18/11*y(i+2) + 9/11*y(i+1) - 2/11*y(i) = 6/11 * h* f(y(i+3))

if method == 3
    disp('Starting BDF3')
     
    for h=hs1

        T = 0.1;
        ts = h:h:T;
        N2 = length(ts);
        y = zeros(N2,N);
        y(1,:) = ones(1,N);

        for i=2:3
            k1 = f(y(i-1,:));
            k2 = f(y(i-1,:) + h/2 .* k1');
            k3 = f(y(i-1,:) + h/2 .* k2');
            k4 = f(y(i-1,:) + h .* k3');
                
            y(i,:) = y(i-1,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4)';
        end

        tic
        A2 = eye(N) + 6/11 * h * A;
        
        fprintf('ConjGrad , h=%d, Nsteps = %d \n',h,N2)
        tol = h^3;
        for i = 1:N2-3
            b2 = 18/11 * y(i+2,:)' - 9/11 * y(i+1,:)' + 2/11 * y(i,:)';
            y(i+3,:) = conjgrad(A2, b2, y(i+3,:)', tol);
        end
        error = y(end,:)'-y_exact;
        fprintf('ConjGrad, Nsteps = %d, error = %d \n', N2, norm(error, inf))
        
        toc
    end
end


%% Conjugate Gradient Method
function x = conjgrad(A, b, x, tol)
    r = b - A * x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        z = A * p;
        alpha = rsold / (z' * p);
        x = x + alpha * p;
        r = r - alpha * z;
        rsnew = r' * r;
        if sqrt(rsnew) < tol
              break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end
