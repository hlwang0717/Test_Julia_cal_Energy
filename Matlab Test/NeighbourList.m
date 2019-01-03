function [ nlist ] = NeighbourList( Config, rcutfact )

Num = Config.Num;
Length = Config.Length; Height = Config.Height;
Rad = Config.Rad;
R_m = rcutfact*max(Rad);

x = [Config.PositX; Config.PositY];
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
    
    % Secondly, calculate the distance between particles.
    % Distances are stored in a matrix Dist, and Dist(i,j) is the distance
    % between particle i and particle j.
        dX = bsxfun(@minus, X_temp, X_temp');
        dY = bsxfun(@minus, Y_temp, Y_temp');
        
        dX = mod(dX+Length/2, Length) - Length/2;
        dY = mod(dY+Height/2, Height) - Height/2;
        
        dx = dX + Strain*dY;
        dy = dY;
        
        Dist = sqrt(dx.^2 + dy.^2);
        Dist(logical(eye(Num))) = Length + Height; %Distance between particle 
    end                                            %i & i is infinite
end

Dist(Dist>R_m) = 0;

[Pair(:,2), Pair(:,1), r] = find(Dist);

% [Pair, id_temp] = sortrows(sort(Pair,2));
% r = r(id_temp);
% Pair = Pair(1:2:end,:);

nlist.i = Pair(:, 1);
nlist.j = Pair(:, 2);
nlist.r = r;

end

