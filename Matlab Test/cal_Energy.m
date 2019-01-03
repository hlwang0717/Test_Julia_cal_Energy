function [ V, dV ] = cal_Energy( x, Config, Pair, NPair )
%CAL_ENERGY1 Summary of this function goes here
%   Detailed explanation goes here

Num = Config.Num;
E0 = Config.E0; Alpha = Config.Alpha;
    
Length = Config.Length; Height = Config.Height;
Rad = Config.Rad;
    
V  = 0;
dV = zeros(2*Num,1);

if Config.BC == 0
    Strain = Config.Strain;
    
    if Config.StrainType == 1
    % Firstly, exert the periodic BC
        x_temp = x(1:Num); y_temp = x(Num+1:end);
        
        X_temp = x_temp - Strain*y_temp;
        Y_temp = y_temp;
        
        X_temp = mod(X_temp, Length);
        Y_temp = mod(Y_temp, Height);
        
        x_temp = X_temp + Strain*Y_temp;
        y_temp = Y_temp;
        
        x = [x_temp; y_temp];
    end

% % %     for k = 1:NPair
% % %         i = Pair(k,1); Xi = X_temp(i); Yi = Y_temp(i); 
% % %         j = Pair(k,2); Xj = X_temp(j); Yj = Y_temp(j);
% % %         
% % %         dX = mod(Xi - Xj + Length/2, Length) - Length/2;
% % %         dY = mod(Yi - Yj + Height/2, Height) - Height/2;
% % %         
% % %         dx = dX + Strain*dY;
% % %         dy = dY;
% % %         
% % %         dij = sqrt(dx^2+dy^2);
% % %         rij = Rad(i) + Rad(j);
% % %         
% % %         temp = rij - dij;
% % %         temp = (abs(temp) + temp)/2;
% % %         
% % %         V = V + E0*temp^Alpha;
% % %         const = -Alpha*E0*temp^(Alpha - 1)/dij;
% % %         
% % %         dV(i) = dV(i) + const*dx; dV(i+Num) = dV(i+Num) + const*dy;
% % %         dV(j) = dV(j) - const*dx; dV(j+Num) = dV(j+Num) - const*dy;
% % %     end
    
    i = Pair(:, 1); j = Pair(:, 2);
    
    Xi = X_temp(i); Yi = Y_temp(i);
    Xj = X_temp(j); Yj = Y_temp(j);
    
    dX = mod(Xi - Xj + Length/2, Length) - Length/2;
    dY = mod(Yi - Yj + Height/2, Height) - Height/2;
    
    dx = dX + Strain*dY;
    dy = dY;
    
    dij = sqrt(dx.^2 + dy.^2);
    rij = Rad(i) + Rad(j);
    
    temp = rij - dij;
    temp = (abs(temp)+temp)/2;

%     rij = Rad(nlist.i) + Rad(nlist.j);
%     temp = rij - nlist.r;
%     temp = (abs(temp)+temp)/2;
    V = V + E0*sum(temp.^Alpha);
%     V = gather(V); % didi-gpu

%     const = -Alpha*E0*temp.^(Alpha-1)./dij;
%     
%     dV(1:Num) = dV(1:Num) + accumarray(i, const.*dx, [Num 1]);
%     dV(1+Num:end) = dV(1+Num:end) + accumarray(i, const.*dy, [Num 1]);
%     dV(1:Num) = dV(1:Num) - accumarray(j, const.*dx, [Num 1]);
%     dV(1+Num:end) = dV(1+Num:end) - accumarray(j, const.*dy, [Num 1]);
%    
    % gpu-array
%     dV = gather(dV);
end

%end

