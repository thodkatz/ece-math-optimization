%%% MINIMIZATION USING DERIVATIVES
clear
clc

%% PLOT OBJECTIVE FUNCTION

f = @(x,y) x^3 * exp(-x^2-y^4);
figure;fsurf(f)
syms x y
fx = matlabFunction(diff(sym(f),x));
fy = matlabFunction(diff(sym(f),y));

initialConditions = [0,0 ; -1,1 ; 1,1];

%% STEEPEST DECENT
% stable learning rate
[fig,k] = steepestDecent(initialConditions(2,:),f,fx,fy,'stable',0.5);

% learning rate that minimize the f(x_k + g_k * d_k)
[fig,k] = steepestDecent(initialConditions(2,:),f,fx,fy,'min',0.5);

% Armijo
[fig,k] = steepestDecent(initialConditions(2,:),f,fx,fy,'armijo',0.5);

%% NEWTON


%% Levenberg Marquardt


% FUNCTIONS
function [fig,k] = steepestDecent(x0,f,fx,fy,gammaFlag,gammaStable)
    fig = figure;fcontour(f,'Fill','on'); colorbar
    hold on

    epsilon = 1e-6;
    maxIterations = 1e3;
    xk = x0(1);
    yk = x0(2);
    plot(xk,yk,'h','color','red','MarkerSize',8)
    %text(xk,yk,"x0",'Color','red','FontSize',14)
    xkNext = @(gamma,xk,yk) xk - gamma * fx(xk,yk);
    ykNext = @(gamma,xk,yk) yk - gamma * fy(xk,yk);
    k = 1;
    while norm(fx(xk,yk),fy(xk,yk)) > epsilon && k<maxIterations
        if gammaFlag == "stable"
            gamma = gammaStable;
        elseif gammaFlag == "min"
            fMin = @(gamma) f(xkNext(gamma,xk,yk),ykNext(gamma,xk,yk));
            gamma = fminsearch(fMin,0);
        elseif gammaFlag == "armijo"
            gamma = armijo(f,fx,fy,xk,yk);
        else 
            disp("Invalid gamma flag")
            exit
        end
        %fprintf("%0.2f,%0.2f\n",xk,yk)
        xk = xkNext(gamma,xk,yk);
        yk = ykNext(gamma,xk,yk);
        k  = k + 1;
        plot(xk,yk,'--o')
    end
    fprintf("Iters: %d/%d\n",k,maxIterations)

    function [gamma] = armijo(f,fx,fy,xk,yk)
        alpha = 1e-1;
        beta  = 0.1;
        s  = 5;
        mk = 0;
        gamma  = s*(beta^mk);
        xkNext = @(gamma,xk,yk) xk - gamma * fx(xk,yk);
        ykNext = @(gamma,xk,yk) yk - gamma * fy(xk,yk);
        xkN = xkNext(gamma,xk,yk);
        ykN = ykNext(gamma,xk,yk);
        D = -fx(xk,yk)^2 -fy(xk,yk)^2;
        while (f(xk,yk) - f(xkN,ykN)) < -alpha*gamma*D
            mk = mk + 1;
            gamma = s*(beta^mk);
            xkN = xkNext(gamma,xk,yk);
            ykN = ykNext(gamma,xk,yk);
        end
    end
end

