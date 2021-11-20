%% Exercise # 1 -  Numerical methods for ODEs
% Course: Numerical Methods For Differential Equations
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% November 2021

close all
clear all


% load("exact_solution.txt")
% t_exact = exact_solution(:,1);
% y_exact = exact_solution(:,2);

y_exact = load("accurate_solution.txt");


%% Question 4

% y' = -Ay;
%% Constants
nx = 30; %100;
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

%% Stability with RK4
% tic
% for i=2:N
%     k1 = f(y(i-1,:));
%     k2 = f(y(i-1,:)+k1'.*h/2);
%     k3 = f(y(i-1,:)+h/2.*k2');
%     k4 = f(y(i-1,:)+h.*k3');
%         
%     y(i,:) = y(i-1,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4)';
% end
% toc

%% Solve with ODE45
% disp('Starting ODE45')
% tic
% tspan = [0 0.1];
% [t,y] = ode45(@(t,y) -A*y, tspan, y(1,:));
% toc
% yf = y(end,:)';
% plot(ts, abs(yf-y_exact))
% return

%% Initial
y = zeros(N,N);
y(1,:) = ones(1,N);

for i=2:3
    k1 = f(y(i-1,:));
    k2 = f(y(i-1,:)+k1'.*h/2);
    k3 = f(y(i-1,:)+h/2.*k2');
    k4 = f(y(i-1,:)+h.*k3');
        
    y(i,:) = y(i-1,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4)';
end

y_init = y;

yi = y(1,:)';
ords=3:5;
hs1=10.^-ords;

%% Solve with CN
disp('Starting CN')

tic
res = zeros(N,3);
iola=1;
for h=hs1
    fprintf('h=%d\n',h)
    A2 = (eye(N)+A*h/2);
    b = (eye(N)-A*h/2);
    % Matlab \ 
    for i=2:N
        fprintf('i=%d\n',i)
        b2 = b*yi;
        yi = A2\b2;
    end
    % Conjugate Gradient:
%     tol = h^3;
%     for i=2:N
%         fprintf('i=%d\n',i)
%         b2 = b*yi;
%         yi = conjgrad(A2,b2,yi,tol);
%     end
    res(:,iola) = yi;
    iola = iola + 1;    
end
toc


%% Solve with BDF3
% y(i+3) - 18/11*y(i+2) + 9/11*y(i+1) - 2/11*y(i) = 6/11 * h* f(y(i+3))

% disp('Starting BDF3')
% tic 
% 
% res = zeros(N,3);
% iola=1;   
% for h=hs1
%     fprintf('h=%d\n',h)
%     A2 = eye(N)+6*h/11*A;
%     y = y_init;
%     % Matlab \
% %     for i=1:N-3
% %         fprintf('i=%d\n',i)
% %         b2 = 18/11*y(i+2,:)' - 9/11*y(i+1,:)' + 2/11*y(i,:)';
% %         y(i+3,:) = A2\b2;
% %     end
%     % Conjugate Gradient:
%     tol = h^3;
%     for i=1:N-3
%         fprintf('i=%d\n',i)
%         b2 = 18/11*y(i+2,:)' - 9/11*y(i+1,:)' + 2/11*y(i,:)';
%         y(i+3,:) = conjgrad(A2,b2,y(i+3,:)',tol);
%     end
%     res(:,iola) = y(end,:);
%     iola = iola + 1;
% end
% toc

%% Plots
figure(1)
yf = res(:,1);
plot(ts, abs(yf-y_exact))
figure(2)
yf = res(:,2);
plot(ts, abs(yf-y_exact))
figure(3)
yf = res(:,3);
plot(ts, abs(yf-y_exact))


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


