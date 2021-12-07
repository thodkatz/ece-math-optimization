%%% MINIMIZATION USING DERIVATIVES
clear
clc

%% PLOT OBJECTIVE FUNCTION

f = @(x,y) x^3 * exp(-x^2-y^4);
%figure;fsurf(f)
%exportgraphics(gcf,'surf.pdf','ContentType','Vector')
%figure;fcontour(f,'Fill','On'); colorbar
%exportgraphics(gcf,'contour.pdf','ContentType','Vector')

initialConditions = [0,0 ; -1,-1 ; 1,1];
epsilon = 1e-6;
maxIterations = 1e3;

%% STEEPEST DECENT
[fig,k] = steepestDescent(maxIterations,epsilon,initialConditions(3,:),f,'stable',0.5);
exportgraphics(gcf,'steepestDescent_(1,1)_stable.pdf','ContentType','Vector')
[fig,k] = steepestDescent(maxIterations,epsilon,initialConditions(3,:),f,'min',0.5);
exportgraphics(gcf,'steepestDescent_(1,1)_min.pdf','ContentType','Vector')
[fig,k] = steepestDescent(maxIterations,epsilon,initialConditions(3,:),f,'armijo',0.5);
exportgraphics(gcf,'steepestDescent_(1,1)_armijo.pdf','ContentType','Vector')

%% NEWTON
[fig,k] = newton(maxIterations,epsilon,initialConditions(3,:),f,'stable',0.1);
exportgraphics(gcf,'newton_(1,1)_stable.pdf','ContentType','Vector')
[fig,k] = newton(maxIterations,epsilon,initialConditions(3,:),f,'min',0.1);
exportgraphics(gcf,'newton_(1,1)_min.pdf','ContentType','Vector')
[fig,k] = newton(maxIterations,epsilon,initialConditions(3,:),f,'armijo',0.1);
exportgraphics(gcf,'newton_(1,1)_armijo.pdf','ContentType','Vector')

%% Levenberg Marquardt
[fig,k] = levenMarq(maxIterations,epsilon,initialConditions(3,:),f,'stable',0.1);
exportgraphics(gcf,'leven_(1,1)_stable.pdf','ContentType','Vector')
[fig,k] = levenMarq(maxIterations,epsilon,initialConditions(3,:),f,'min',0.1);
exportgraphics(gcf,'leven_(1,1)_min.pdf','ContentType','Vector')
[fig,k] = levenMarq(maxIterations,epsilon,initialConditions(3,:),f,'armijo',0.1);
exportgraphics(gcf,'leven_(1,1)_armijo.pdf','ContentType','Vector')


% FUNCTIONS
function [fig,k] = steepestDescent(maxIterations,epsilon,x0,f,gammaFlag,gammaStable)
    fig = figure;fcontour(f,'Fill','on'); colorbar
    hold on

    x = sym('x');
    y = sym('y');
    fx = matlabFunction(diff(sym(f),x));
    fy = matlabFunction(diff(sym(f),y));

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
            gamma = fminbnd(fMin,0,5);
        elseif gammaFlag == "armijo"
            initStep = 2;
            [gamma,mk] = armijo(initStep,f,fx,fy,@(x,y)-fx(x,y),@(x,y)-fy(x,y),xk,yk,xkNext,ykNext);
        else 
            disp("Invalid gamma flag")
            exit
        end
        %fprintf("%0.2f,%0.2f,gamma=%d\n",xk,yk,gamma)
        xk = xkNext(gamma,xk,yk);
        yk = ykNext(gamma,xk,yk);
        k  = k + 1;
        plot(xk,yk,'--o')
    end
    fprintf("Iters: %d/%d\n",k,maxIterations)
end

function [fig,k] = newton(maxIterations,epsilon,x0,f,gammaFlag,gammaStable)
    fig = figure;fcontour(f,'Fill','on'); colorbar
    hold on

    x = sym('x');
    y = sym('y');
    fx  = diff(sym(f),x);
    fy  = diff(sym(f),y);
    fxx = diff(sym(fx),x);
    fxy = diff(sym(fx),y);
    fyx = diff(sym(fy),x);
    fyy = diff(sym(fy),y);
    J = [fxx fxy ; fyx fyy];
    dk  = -inv(J) * [fx ; fy];
    dkX = matlabFunction(dk(1));
    dkY = matlabFunction(dk(2));
    fx  = matlabFunction(sym(fx));
    fy  = matlabFunction(sym(fy));
    xkNext = @(gamma,xk,yk) xk + gamma * dkX(xk,yk);
    ykNext = @(gamma,xk,yk) yk + gamma * dkY(xk,yk);


    xk = x0(1);
    yk = x0(2);
    plot(xk,yk,'h','color','red','MarkerSize',8)
    %text(xk,yk,"x0",'Color','red','FontSize',14)
    k = 1;
    while norm(fx(xk,yk),fy(xk,yk)) > epsilon && k<maxIterations
        if gammaFlag == "stable"
            gamma = gammaStable;
        elseif gammaFlag == "min"
            fMin = @(gamma) f(xkNext(gamma,xk,yk),ykNext(gamma,xk,yk));
            gamma = fminbnd(fMin,0,5);
        elseif gammaFlag == "armijo"
            initStep = 2;
            [gamma,mk] = armijo(initStep,f,fx,fy,dkX,dkY,xk,yk,xkNext,ykNext);
        else 
            disp("Invalid gamma flag")
            exit
        end
        %fprintf("%0.2f,%0.2f,gamma=%d\n",xk,yk,gamma)
        xk = xkNext(gamma,xk,yk);
        yk = ykNext(gamma,xk,yk);
        k  = k + 1;
        plot(xk,yk,'--o')
    end
    fprintf("Iters: %d/%d\n",k,maxIterations)
end

function [fig,k] = levenMarq(maxIterations,epsilon,x0,f,gammaFlag,gammaStable)
    fig = figure;fcontour(f,'Fill','on'); colorbar
    hold on

    x = sym('x');
    y = sym('y');
    fx  = diff(sym(f),x);
    fy  = diff(sym(f),y);
    fxx = diff(sym(fx),x);
    fxy = diff(sym(fx),y);
    fyx = diff(sym(fy),x);
    fyy = diff(sym(fy),y);
    J = [fxx fxy ; fyx fyy];
    J = matlabFunction(sym(J));
    fx = matlabFunction(diff(sym(f),x));
    fy = matlabFunction(diff(sym(f),y));

    xk = x0(1);
    yk = x0(2);
    plot(xk,yk,'h','color','red','MarkerSize',8)
    %text(xk,yk,"x0",'Color','red','FontSize',14)
    k = 1;
    while norm(fx(xk,yk),fy(xk,yk)) > epsilon && k<maxIterations
        hessian = J(xk,yk);
        eigen = eig(hessian);
        m = abs(max(eigen))+0.1;
        d = - inv(hessian+m*eye(2)) * [fx(xk,yk) ; fy(xk,yk)];

        if gammaFlag == "stable"
            gamma = gammaStable;
        elseif gammaFlag == "min"
            fMin = @(gamma) f(xk+gamma*d(1),yk+gamma*d(2));
            gamma = fminbnd(fMin,0,5);
        elseif gammaFlag == "armijo"
            initStep = 2;
            [gamma,mk] = armijo(initStep,f,fx,fy,J,xk,yk);
        else 
            disp("Invalid gamma flag")
            exit
        end

        %fprintf("%0.2f,%0.2f,gamma=%d\n",xk,yk,gamma)
        xk = xk + gamma * d(1);
        yk = yk + gamma * d(2);
        k  = k + 1;
        plot(xk,yk,'--o')
    end
    fprintf("Iters: %d/%d\n",k,maxIterations)

    % TODO: Refactoring. Concistency among armijo functions
    function [gamma,mk] = armijo(s,f,fx,fy,J,xk,yk)
        alpha = 1e-1;
        beta  = 0.1;
        mk = 0;
        gamma  = s*(beta^mk);

        hessian = J(xk,yk);
        eigen = eig(hessian);
        m = abs(max(eigen))+0.1;
        d = - inv(hessian+m*eye(2)) * [fx(xk,yk) ; fy(xk,yk)];
        xkN = xk + gamma * d(1);
        ykN = yk + gamma * d(2);

        D = fx(xk,yk)*d(1) + fy(xk,yk)*d(2);
        while (f(xk,yk) - f(xkN,ykN)) < -alpha*gamma*D
            mk = mk + 1;
            gamma = s*(beta^mk);

            hessian = J(xk,yk);
            eigen = eig(hessian);
            m = abs(max(eigen))+0.1;
            d = - inv(hessian+m*eye(2)) * [fx(xk,yk) ; fy(xk,yk)];
            xkN = xk + gamma * d(1);
            ykN = yk + gamma * d(2);
         end
    end
end

function [gamma,mk] = armijo(s,f,fx,fy,dkX,dkY,xk,yk,xkNext,ykNext)
    alpha = 1e-1;
    beta  = 0.1;
    mk = 0;
    gamma  = s*(beta^mk);
    xkN = xkNext(gamma,xk,yk);
    ykN = ykNext(gamma,xk,yk);
    D = fx(xk,yk)*dkX(xk,yk) + fy(xk,yk)*dkY(xk,yk);
    while (f(xk,yk) - f(xkN,ykN)) < -alpha*gamma*D
        mk = mk + 1;
        gamma = s*(beta^mk);
        xkN = xkNext(gamma,xk,yk);
        ykN = ykNext(gamma,xk,yk);
    end
end
