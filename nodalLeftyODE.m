
%% [nodal lefty leftyInhibitor]
% f1 = @(uu1,uu2,uu3,uu4) (param.K11_1).*(uu1.^param.n11)./((param.K11_2^param.n11 + uu1.^param.n11).*(1+param.K21.*uu2)) ...
%     - (param.kd1).*uu1 + param.K41_1.*(uu4.^param.n41)./(param.K41_2.^param.n41+ uu4.^param.n41);
% f2 = @(uu1,uu2,uu3,uu4) (param.K12_1).*(uu1.^param.n12)./(param.K12_2^param.n12 + uu1.^param.n12)...
%     - (param.kd2).*uu2 - (param.K32).*uu3;
% f3 = @(uu1,uu2,uu3,uu4) (param.K13_1).*(uu1.^param.n13)./(param.K13_2^param.n13 + uu1.^param.n13) - (param.kd3).*uu3 ;

%% Parameters (p)
% K11_1 (1)
% K11_2 (2)
% n11   (3)
% K21   (4)
% Kd1   (5)
% 
% K12_1 (6)
% K12_2 (7)
% n12   (8)
% K32   (9)
% Kd2   (10)
% 
% K13_1 (11)
% K13_2 (12)
% n13   (13)
% Kd3   (14)


%%

function dYdt = nodalLeftyODE(t,Y,p)
dYdt = zeros(3,1);
dYdt(1) = p(1)*(Y(1)^p(3))./((p(2)^p(3) + Y(1)^p(3)).*(1+p(4)*Y(2))) - p(5)*Y(1); 
dYdt(2) = p(6)*(Y(1)^p(8))./(p(7)^p(8) + Y(1)^p(8)) - p(9)*Y(3) - p(10)*Y(2);
dYdt(3) = p(11)*(Y(1)^p(13))./(p(12)^p(13) + Y(1)^p(13)) - p(14)*Y(3);
