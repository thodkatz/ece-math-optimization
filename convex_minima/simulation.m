%% MINIMIZATION PROBLEM FOR CONVEX FUNCTIONS SINGLE VARIABLE
clc
clear all

% FUNCTIONS
f1 = @(x) (x-3)^2 + (sin(x+3))^2;
f2 = @(x) (x-1) * cos(0.5 * x) + x^2;
f3 = @(x) (x+2)^2 + exp(x-2) * sin(x+3);
numFuncs = 3;

intervalStart = [-4 4];

syms x;
figure
hold on
fplot(f1(x),intervalStart);
fplot(f2(x),intervalStart);
fplot(f3(x),intervalStart);

%% BISECTION METHOD
% STABLE ACCURACY (lambda) VARIABLE STEP (epsilon)
intervalAccuracy = 0.01;
intervalStep = [0.0001 0.0004 0.0006 0.0008 0.001 0.002 0.003 0.004 0.0042 0.0046 0.0048];
ends = cell(numFuncs,numel(intervalStep));
countFunctionCalls = cell(numFuncs,numel(intervalStep));

funcs=cell(numFuncs,1);
funcs{1} = f1; funcs{2} = f2; funcs{3} = f3;
for idFunc=1:numel(funcs)
    for i=1:numel(intervalStep)
        [endPoints,count] = bisection(intervalStart,intervalStep(i),intervalAccuracy,funcs{idFunc});
        ends(idFunc,i) = {endPoints};
        countFunctionCalls(idFunc,i) = {count};
    end
end

%% BISECTION METHOD
% STABLE STEP (epsilon) VARIABLE ACCURACY (lambda)
intervalStep = 0.001;
intervalAccuracy = [0.0025 0.004 0.006 0.008 0.01 0.012 0.014 0.016 0.018 0.02 0.025];
ends = cell(numFuncs,numel(intervalAccuracy));
countFunctionCalls = cell(numFuncs,numel(intervalAccuracy));

funcs=cell(numFuncs,1);
funcs{1} = f1; funcs{2} = f2; funcs{3} = f3;
for idFunc=1:numel(funcs)
    for i=1:numel(intervalAccuracy)
        [endPoints,count] = bisectionMethod(intervalStart,intervalStep,intervalAccuracy(i),funcs{idFunc});
        ends(idFunc,i) = {endPoints};
        countFunctionCalls(idFunc,i) = {count};
    end
end

%% GOLDEN SECTOR METHOD
intervalAccuracy = [0.0025 0.004 0.006 0.008 0.01 0.012 0.014 0.016 0.018 0.02 0.025];
ends = cell(numFuncs,numel(intervalAccuracy));
countFunctionCalls = cell(numFuncs,numel(intervalAccuracy));

funcs=cell(numFuncs,1);
funcs{1} = f1; funcs{2} = f2; funcs{3} = f3;
for idFunc=1:numel(funcs)
    for i=1:numel(intervalAccuracy)
        [endPoints,count] = goldenSector(intervalStart,intervalAccuracy(i),funcs{idFunc});
        ends(idFunc,i) = {endPoints};
        countFunctionCalls(idFunc,i) = {count};
    end
end

%% FIBONACCI METHOD
intervalStep = 0.0001;
intervalAccuracy = [0.0025 0.004 0.006 0.008 0.01 0.012 0.014 0.016 0.018 0.02 0.025];
ends = cell(numFuncs,numel(intervalAccuracy));
countFunctionCalls = cell(3,numel(intervalAccuracy));

funcs=cell(numFuncs,1);
funcs{1} = f1; funcs{2} = f2; funcs{3} = f3;
for idFunc=1:numel(funcs)
    for i=1:numel(intervalAccuracy)
        [endPoints,count] = fibonacci(intervalStart,intervalStep,intervalAccuracy(i),funcs{idFunc});
        ends(idFunc,i) = {endPoints};
        countFunctionCalls(idFunc,i) = {count};
    end
end

%% PLOTS
% varying epsilon
figure
plot(intervalStep,cell2mat(countFunctionCalls(1,:)),'--o')
xlabel('$\epsilon$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_1(x)$','fontsize',14,'interpreter','latex')

figure
plot(intervalStep,cell2mat(countFunctionCalls(2,:)),'--o')
xlabel('$\epsilon$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_2(x)$','fontsize',14,'interpreter','latex')

figure
plot(intervalStep,cell2mat(countFunctionCalls(3,:)),'--o')
xlabel('$\epsilon$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_3(x)$','fontsize',14,'interpreter','latex')

%% varying end points with stable lambda
figure
hold on
intervalStepIndex = 11;
for i = 1:3
pairs = cell2mat(ends(i,intervalStepIndex));
[k,] = size(pairs);
plot(1:k,pairs(:,1),'--o')
plot(1:k,pairs(:,2),'--o')
xlabel('$\epsilon$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_1(x)$','fontsize',14,'interpreter','latex')
end

%% varying lambda
figure
plot(intervalAccuracy,cell2mat(countFunctionCalls(1,:)),'--o')
xlabel('$\lambda$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_1(x)$','fontsize',14,'interpreter','latex')

figure
plot(intervalAccuracy,cell2mat(countFunctionCalls(2,:)),'--o')
xlabel('$\lambda$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_2(x)$','fontsize',14,'interpreter','latex')

figure
plot(intervalAccuracy,cell2mat(countFunctionCalls(3,:)),'--o')
xlabel('$\lambda$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_3(x)$','fontsize',14,'interpreter','latex')

%% varying end points with stable epsilon
figure
hold on
intervalAccuracyIndex = 5;
for i = 1:3
pairs = cell2mat(ends(i,intervalAccuracyIndex));
[k,] = size(pairs);
plot(1:k,pairs(:,1),'--o')
plot(1:k,pairs(:,2),'--o')
xlabel('$\epsilon$','fontsize',14,'interpreter','latex')
ylabel('Function calls $f_1(x)$','fontsize',14,'interpreter','latex')
end

function [ends,countFunctionCalls] = bisection(intervalStart,intervalStep,intervalAccuracy,f)
    a(1) = intervalStart(1);
    b(1) = intervalStart(2);
    countFunctionCalls = 0;
    i = 1;
    while b(i)-a(i) >= intervalAccuracy
        midPoint = (a(i)+b(i))/2;
        x1 = midPoint - intervalStep;
        x2 = midPoint + intervalStep;
        if f(x1) < f(x2)
            a(i+1) = a(i);
            b(i+1) = x2;
        else
            a(i+1) = x1;
            b(i+1) = b(i);
        end
        countFunctionCalls = countFunctionCalls+2;
        i = i+1;
    end
    
    ends(:,1) = a;
    ends(:,2) = b;
end

function [ends,countFunctionCalls] = goldenSector(intervalStart,intervalAccuracy,f)
    % init
    gamma = 0.618;
    a(1) = intervalStart(1);
    b(1) = intervalStart(2);
    nextX1 = @(a,b,gamma) a + (1-gamma)*(b-a);
    nextX2 = @(a,b,gamma) a + gamma*(b-a);
    x1 = nextX1(a(1),b(1),gamma);
    x2 = nextX2(a(1),b(1),gamma);
    y1 = f(x1);
    y2 = f(x2);
    
    i = 1;
    countFunctionCalls = 2;  
    while b(i)-a(i) >= intervalAccuracy
        if y1 < y2
            a(i+1) = a(i);
            b(i+1) = x2;
            x2 = x1;
            x1 = nextX1(a(i+1),b(i+1),gamma);
            y2 = y1;
            y1 = f(x1);
        else
            a(i+1) = x1;
            b(i+1) = b(i);
            x1 = x2;
            x2 = nextX2(a(i+1),b(i+1),gamma);
            y1 = y2;
            y2 = f(x2);
        end
        countFunctionCalls = countFunctionCalls+1;
        i = i+1;
    end
    
    ends(:,1) = a;
    ends(:,2) = b;
end

function [ends,countFunctionCalls] = fibonacci(intervalStart,intervalStep,intervalAccuracy,f)
    % init
    gamma = 0.618;
    a(1) = intervalStart(1);
    b(1) = intervalStart(2);
    nextX1 = @(a,b,gamma) a + (1-gamma)*(b-a);
    nextX2 = @(a,b,gamma) a + gamma*(b-a);
    x1 = nextX1(a(1),b(1),gamma);
    x2 = nextX2(a(1),b(1),gamma);
    y1 = f(x1);
    y2 = f(x2);
    
    i = 1;
    countFunctionCalls = 2;  
    while b(i)-a(i) >= intervalAccuracy
        if y1 < y2
            a(i+1) = a(i);
            b(i+1) = x2;
            x2 = x1;
            x1 = nextX1(a(i+1),b(i+1),gamma);
            y2 = y1;
            y1 = f(x1);
        else
            a(i+1) = x1;
            b(i+1) = b(i);
            x1 = x2;
            x2 = nextX2(a(i+1),b(i+1),gamma);
            y1 = y2;
            y2 = f(x2);
        end
        countFunctionCalls = countFunctionCalls+1;
        i = i+1;
    end
    
    ends(:,1) = a;
    ends(:,2) = b;
end