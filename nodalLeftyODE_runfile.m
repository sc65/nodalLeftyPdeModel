


K11_1 = 0.11;    %(1)
K11_2 = 1;
n11   = 5;
K21   = 0.1;
Kd1   = 0.001;  %(5)

K12_1 = 0.02;
K12_2 = 0.5;
n12   = 5;      
K32   = 0.03;    
Kd2   = 0.001;   %(10)

K13_1 = 0.01;    
K13_2 = 12;     
n13   = 5;       
Kd3   = 0.0002;  %(14)


p = [K11_1, K11_2, n11, K21, Kd1,...
    K12_1, K12_2, n12, K32, Kd2,...
    K13_1, K13_2, n13, Kd3];
    
%%

tspan = [0 350];
Y0 = [32 2 0.001];
[t,y] = ode45(@(t,y) nodalLeftyODE(t,y,p), tspan, Y0);
%[t,y] = ode15s(@(t,y) nodalLeftyODE(t,y,p), tspan, Y0);

figure; 
for ii = 1:3
subplot(1,3,ii);
plot(t, y(:,ii), 'k-');
end
%%

%% ---------- sensitivity analyses
[t,y] = ode45(@(t,y) nodalLeftyODE(t,y,p), tspan, Y0);
solution1 = y;

pCost = zeros(numel(p),3);
p1 = p;
figure; hold on;
for ii = 1:numel(p)
    p1(ii) = p(ii)*2;
    [t,y] = ode45(@(t,y) nodalLeftyODE(t,y,p1), tspan, Y0);
    
    subplot(3,5,ii);
    plot(t, y(:,2), 'k-');
    title(num2str(ii));
    %pCost(ii,:) = sum((solution1-y).^2,1);
    p1 = p;
end

    
   
%% ------------- finding steady states
fun = @(y) findSteadyStates(y, p);
Y0 = [200 ; 2; 0.1];
options = optimoptions('fsolve','PlotFcn',@optimplotfirstorderopt, 'MaxFunctionEvaluations', 1000);
ss = fsolve(fun, Y0, options);

function F = findSteadyStates(Y, p)
% solves the system of equations F(Y) = 0;
F(1) = p(1)*(Y(1)^p(3))./((p(2)^p(3) + Y(1)^p(3)).*(1+p(4)*Y(2))) - p(5)*Y(1); 
F(2) = p(6)*(Y(1)^p(8))./(p(7)^p(8) + Y(1)^p(8)) - p(9)*Y(3) - p(10)*Y(2);
F(3) = p(11)*(Y(1)^p(13))./(p(12)^p(13) + Y(1)^p(13)) - p(14)*Y(3);
end
