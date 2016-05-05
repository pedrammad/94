% May/05/2016
%UCLA MAE 94
% COPY RIGHT, PEDRAM MADANCHIAN

close all; clear all; clc


x_constraint = 0.2;              % [m]  length
z_constraint = 0.1;              % [m]  width
y_constraint = 0.1;              % [m]  height

k_spring = 17.38;                % [lbf/in] Spring constant
g = 32.2*12;                        % [in/s^2]
rho_ball = 0.283;                % [lbm/in^3] 
r_ball = 3/8;                    % [in]
vol_ball = 4/3 * pi *r_ball^3;
m_ball =rho_ball*vol_ball ;      % [kg]

x = 0.1:0.001:x_constraint;      % [m]
x = x * 39.3701;                 % [in]
y = 0.05:0.001:y_constraint;      % [m]
y = y * 39.3701;                 % [in]

totalCanons = length(y)*length(x);
canon = zeros(totalCanons, 2);
iCanon = 0;

for i = 1:length(x)
    for j = 1:length(y)
        iCanon = iCanon + 1;
        canon (iCanon, 1) = x(i);
        canon (iCanon, 2) = y(j);
        

    end
end

for k = 1:iCanon
   x (k) = canon(k,1);
   y (k) = canon(k,2);
   
   alpha(k) = atan(y(k)./x(k));        %[rad] found based on x and y
   L(k) = x(k)./ cos(alpha(k));        %[in] Canon Length
   v_top(k) = sqrt (2*g*L(k).*sin(alpha(k)) + k_spring/m_ball*(L(k)-2).^2); %[in/s]
   
   %[in] Range Equation
   d(k) = v_top(k).*cos(alpha(k))./g *(v_top(k).*sin(alpha(k)) + ...
          sqrt((v_top(k).*sin(alpha(k)))^2 + 2*g*L(k).*sin(alpha(k))));
   alpha (k) = alpha(k) * 180/pi;
end

[maxValue,maxIndex] = max(d);

fprintf (['Maximum Range of ' num2str(d(maxIndex)/39.3701) ' inches occurs @ Case ' ...
                  num2str(maxIndex) ' When: \n Length  =\t' ...
                  num2str(x(maxIndex)/39.3701) ' meters & \n Height =\t' ...
                  num2str(y(maxIndex)/39.3701) ' meters & \n Alpha =\t' ...
                  num2str(alpha(maxIndex)) ' degrees\n\n']);
              

