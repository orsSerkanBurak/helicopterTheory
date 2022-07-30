clc; clear; close all;
g=9.81; % Gravitional Acceleation [m/sec^2]
rho= 1.225; % Density of Air [kg/m^3]
GTOW=4000; % Aerospatiale AS 365N Gross Take Off Mass [kg];
d = 11.93; % Diameter of Main Rotor [m]
A=pi*(d/2)^2; % Actuator Disk Area [m^2];
T = GTOW*g; % Thrust [N]
V_inf = 0:0.01:90; % Forward Flight Speed [m/s^2], must be less than 90 m/sec^2
omega = 300*(1/60)*(2*pi); % Main Rotor Angular Speed [rad/sec]
Vtip = omega*(d/2); % Tip Speed of Main Rotor
Ct = T/(rho*A*Vtip^2); % Thrust Coefficient in Hover Condition
lambda_h = sqrt(Ct/2); % Initial Inflow Ratio
lambdafp(1) = sqrt(Ct/2);
alphaValues = -2:2:8; % Angle of Attack Values
error = 0.0005; % Permitted Error Percent
N = 1000; % Iteration Number
count1 = 1; % Counting Variable for Fixed Point Iteration Method
count2 = 1; % Counting Variable for Newton-Raphson Method


%% Fixed Point Iteration Method
for i=1:length(alphaValues)
    alpha = alphaValues(i);
    for j = 1:length(V_inf)
        Vinf = V_inf(j);
        Tinf = T/cosd(alpha);
        Ct_inf = Tinf/(rho*A*Vtip^2);
        y = sqrt(Ct_inf/2);
        mufp = (Vinf*cosd(alpha))/Vtip;
        x0 = lambda_h;
        lambda = x0;
        Xold = x0;
        for k = 0:7
            lambda = (mufp*tand(alpha))+(Ct_inf/(2*sqrt(mufp^2+lambda^2)));
            ifp(k+1) = k;
            if abs((lambda-Xold)/lambda) < error
                break
            end
            Xold = lambda;
        end
        a = sqrt(Ct/2);
        lambdaa(count1) = lambda/a;
        muu(count1) = mufp/a;
        count1 = count1+1;
    end
    figure(1)
    plot(muu,lambdaa);
    hold on;
    grid on;
    count1 = 1;
end
title('Fixed Point Iteration Method');
xlabel('\mu/\lambda_{h}');
ylabel('\lambda/\lambda_{h}');
legend('\alpha = -2','\alpha = 0','\alpha = 2','\alpha = 4','\alpha = 6','\alpha = 8','Location','southwest');

%% Newton-Raphson Method
for i=1:length(alphaValues)
    alpha = alphaValues(i);
    for j = 1:length(V_inf)
        Vinf2 = V_inf(j);
        Tinf2 = T/cosd(alpha);
        Ct_inf2 = Tinf2/(rho*A*Vtip^2);
        y12 = sqrt(Ct_inf2/2);
        mu2 = (Vinf2*cosd(alpha))/Vtip;
        f = @(lambda2) lambda2 - mu2*tand(alpha) - Ct_inf2/(2*sqrt(mu2^2+lambda2^2));
        df=@(lambda2) 1+(Ct_inf2/2)*(mu2^2+lambda2^2)^(-3/2)*lambda2;
        for k = 0:N
            f0 = f(lambda_h);
            f0_diff = df(lambda_h);
            y12 = lambda_h - (f0/f0_diff);
            ins(k+1)=k;
            if abs((y12-lambda_h)/y12) < error
                break
            end
            lambda_h = y12;
        end
        a2 = sqrt(Ct/2);
        lambdaa2(count2) = y12/a2;
        muu2(count2) = mu2/a2;
        count2 = count2+1;
    end
    figure(2)
    plot(muu2,lambdaa2);
    hold on;
    grid on;
    count2 = 1;
end
title('Newton Raphson Method');
xlabel('\mu/\lambda_{h}');
ylabel('\lambda/\lambda_{h}');
legend('\alpha = -2','\alpha = 0','\alpha = 2','\alpha = 4','\alpha = 6','\alpha = 8','Location','southwest');

fprintf('Inflow Ratio is converged at %d iteration number with applying Fixed Point Iteration Method\n', ifp(end));
fprintf('Inflow Ratio is converged at %d iteration number with applying Newton-Raphson Method\n', ins(end));