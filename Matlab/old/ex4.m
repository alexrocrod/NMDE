%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all

%% Question 4

% load("exact_solution.txt")
% t_exact = exact_solution(:,1);
% y_exact = exact_solution(:,2);

y_exact = load("accurate_solution.txt");

use_conjgrad=true;

method = 1; % 0 - RK4 stability,  1 - ODE45, 2 - Crank-Nicolson, 3 - BDF3

% y' = -Ay;
%% Constants
nx = 100; %100;

G = numgrid ( 'S' , nx ) ;
A = delsq(G) * ( nx-1)^2 ;
lambda = -eigs(A,1,'lm');

T = 0.1;
h=T/length(A);
ts = h:h:T;

N = length(ts);
y = zeros(N,N);
y(1,:) = ones(1,N);

if nx ~= 100
    disp('Starting exact solution computation')
    tic
    y_exact = expm(-0.1*A)*ones(N,1);
    toc
end
f = @(x) -A*x';

ords = 3:5;
hs1 = 10.^-ords;

%% Stability with RK4
if method == 0
    tic
    for i=2:N
        k1 = f(y(i-1,:));
        k2 = f(y(i-1,:)+k1'.*h/2);
        k3 = f(y(i-1,:)+h/2.*k2');
        k4 = f(y(i-1,:)+h.*k3');
            
        y(i,:) = y(i-1,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4)';
    end
    toc
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
        ola1 = now;
        datetime(ola1,'ConvertFrom','datenum')
        T = 0.1;
        ts = h:h:T;
        N2 = length(ts);
        y = zeros(N2,N);
        y(1,:) = ones(1,N);
        yi = y(1,:)';

        tic
        A2 = eye(N) + A * h/2;
        b = eye(N) - A * h/2;
        if ~use_conjgrad
             % Matlab \ :
            fprintf('Matlab \\ , h=%d, Nsteps = %d \n', h, N2)
            for i=2:N2
                b2 = b * yi;
                yi = A2\b2;
            end
            error = yi-y_exact;
            fprintf('Matlab \\, Nsteps = %d, error = %d \n', N2, norm(error, inf))
        else
            % Conjugate Gradient:
            fprintf('ConjGrad , h=%d, Nsteps = %d \n', h, N2)
            tol = h^3;
            for i=2:N2
                b2 = b * yi;
                yi = conjgrad(A2, b2, yi, tol);
            end
            error = yi-y_exact;
            fprintf('ConjGrad, Nsteps = %d, error = %d \n', N2, norm(error, inf))
        end
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
        if ~use_conjgrad
            fprintf('Matlab \\, h=%d, Nsteps = %d \n',h,N2)
            for i = 1:N2-3
                b2 = 18/11 * y(i+2,:)' - 9/11 * y(i+1,:)' + 2/11 * y(i,:)';
                y(i+3,:) = A2\b2;
            end
            error = y(end,:)'-y_exact;
            fprintf('Matlab \\, Nsteps = %d, error = %d \n', N2, norm(error, inf))
        else
            fprintf('ConjGrad , h=%d, Nsteps = %d \n',h,N2)
            tol = h^3;
            for i = 1:N2-3
                b2 = 18/11 * y(i+2,:)' - 9/11 * y(i+1,:)' + 2/11 * y(i,:)';
                y(i+3,:) = conjgrad(A2, b2, y(i+3,:)', tol);
            end
            error = y(end,:)'-y_exact;
            fprintf('ConjGrad, Nsteps = %d, error = %d \n', N2, norm(error, inf))
        end
        toc
    end
end


%% Conjugate Gradient Method
function x = conjgrad(A, b, x, tol)
    r = b - A * x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < tol
              break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end
