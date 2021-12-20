%%% STEEPEST DESCENT WITH PROJECTION
%clc
clear

syms x y
f(x,y) = 0.5 * x^2 + 2* y^2;
gradf = gradient(f,[x,y]);
gradf = matlabFunction(gradf);
gradfvect = @(x) gradf(x(1),x(2));
f = matlabFunction(f);
fvect = @(x) f(x(1),x(2));

%figure; fsurf(f);
%exportgraphics(gcf,'surf.pdf','ContentType','Vector')

% goal: minimum of f, which is 0
xMin = [0,0];
maxIter = 150;

%% EXERCISE 1
e = 0.01;
x0 = [1;1];
%gammas = [0.05,0.5,2,10];
gammas = [0.05];

for i=1:numel(gammas)
    figure
    xk = steepestDescent(maxIter,e,x0,gradfvect,gammas(i));
    plot(xk(:,1),xk(:,2),'--o')
    xlabel('x1')
    ylabel('x2')
    %exportgraphics(gcf,['surf_' int2str(gammas(i)) '.pdf'],'ContentType','Vector')
end

for i=1:numel(gammas)
    figure
    xk = steepestDescentCondition(maxIter,e,x0,gradfvect,gammas(i));
    plot(xk(:,1),xk(:,2),'--o')
    xlabel('x1')
    ylabel('x2')
    %exportgraphics(gcf,['surf_constraint_' int2str(gammas(i)) '.pdf'],'ContentType','Vector')
end

% %% EXERCISE 2
e = 0.01;
gamma = 0.05;
sk = 8;
x0 = [10;-5];
figure
xk = steepestDescentProjection(maxIter,e,sk,x0,gradfvect,gamma);
myPlot(xk,f,fvect)
%exportgraphics(gcf,'exe2.pdf','ContentType','Vector')

% %% EXERCISE 3
e = 0.02;
gamma = 0.3;
sk = 10;
x0 = [-7;5];
figure
xk = steepestDescentProjection(maxIter,e,sk,x0,gradfvect,gamma);
myPlot(xk,f,fvect)
%exportgraphics(gcf,'exe3.pdf','ContentType','Vector')


% %% EXERCISE 4
e = 0.01;
gamma = 0.1;
sk = 0.5;
x0 = [17;-5];
figure
xk = steepestDescentProjection(maxIter,e,sk,x0,gradfvect,gamma);
myPlot(xk,f,fvect)
exportgraphics(gcf,'exe4.pdf','ContentType','Vector')

function myPlot(xk,f,fvect)
    [k,dims] = size(xk);
    for i=1:k
        fvalues(i) = fvect(xk(i,:));
    end
    subplot(2,1,1)
    plot(1:k,fvalues)
    xlabel('k')
    ylabel('f')
    %subplot(2,2,2)
    %plot(xk(:,1),xk(:,2),'--o')
    subplot(2,1,2)
    hold on
    fcontour(f,[-20 20])
    plot(xk(:,1),xk(:,2),'--o')
    xlabel('x1')
    ylabel('x2')
end

function x = steepestDescent(maxIter,e,x0,gradfvect,gamma)
    k = 1;
    xk = x0;
    x(k,:) = xk; 

    while norm(gradfvect(xk)) > e && k<maxIter
        xk = xk - gamma .* gradfvect(xk);
        k  = k + 1;
        x(k,:) = xk;
    end
    fprintf("Iters: %d/%d\n",k,maxIter)
end

function x = steepestDescentCondition(maxIter,e,x0,gradfvect,gamma)
    k = 1;
    xk = x0;
    x(k,:) = xk;
    while norm(gradfvect(xk)) > e && k<maxIter
        options = optimoptions('fmincon','Display','off');
        constraintPoint = fmincon(@(x)gradfvect(xk)' * ([x(1)-xk(1);x(2)-xk(2)]),[1,1],[],[],[],[],[-15,-20],[15,12],[],options);
        constraintPoint = constraintPoint';
        xk = xk + gamma .* (constraintPoint - xk);
        k = k + 1;
        x(k,:) = xk;
    end
    fprintf("Iters: %d/%d\n",k,maxIter)
end

function x = steepestDescentProjection(maxIter,e,sk,x0,gradfvect,gamma)
    k = 1;
    xk = x0;
    x(k,:) = xk;
    %while ~isCritical(gradfvect,sk,xk) && k<maxIter
    while norm(gradfvect(xk))>e && k<maxIter
        xProjection = project(xk,sk,gradfvect,[-15;15],[-20;12]);
        xk = xk + gamma .* (xProjection - xk);
        k = k + 1;
        x(k,:) = xk;
    end
    fprintf("Iters: %d/%d\n",k,maxIter)

    function out = isCritical(gradfvect,s,x)
        xProjection = x - s * gradfvect(x);
        out = isequal(xProjection,x);
    end

    function xProjection = project(xk,sk,gradfvect,x1bounds,x2bounds)
        x1min = x1bounds(1);
        x1max = x1bounds(2);
        x2min = x2bounds(1);
        x2max = x2bounds(2);
        xProjection = xk - sk * gradfvect(xk);
        if xProjection(1) <= x1min
            xProjection(1) = x1min;
        end
        if xProjection(1) >= x1max
            xProjection(1) = x1max;
        end
        if xProjection(2) <= x2min
            xProjection(2) = x2min;
        end
        if xProjection(2) >= x2max
            xProjection(2) = x2max;
        end
    end
end