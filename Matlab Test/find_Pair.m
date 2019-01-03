function [Pair, contact_number, rattlers] = find_Pair( Config )
% Input:
%   Config : all information of a configuration is stored in structure
%             Config, including granular number, position, radius, boundary
%             condition, etc.

% Output:
%   Pair : a sparse matrix that stored all possible contact pair
%                  information among all particles.
%   rattlers: a structure stores the 
%               rattlers' id (has no overlap withother particles);
%               gap (the minimial distance with other particle);
%               neighbor (particles around the rattlers)

Num = Config.Num;
Length = Config.Length; Height = Config.Height;
Rad = Config.Rad;
R_m = 5*max(Rad);

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
        Dist(logical(eye(Num))) = Length + Height; %Distance between particle i & i is infinite
    end
end

if nargout > 1
    R = bsxfun(@plus, Rad, Rad');
    R_min = min(Rad);
    temp1 = Dist - R;  % temp1(i,j) = Dist(i,j) - (r_i + r_j)
    temp2 = sum(temp1 < 0);
    
    contact_number = temp2';

    rattlers.id = find(temp2 <= 2);

    rattlers.gap = min(temp1(:,rattlers.id));
    
    rattlers.num = length(rattlers.id);

    for i = 1:rattlers.num
        rattlers.neighbor{i} = find(temp1(rattlers.id(i),:)<0.5*R_min);
    end
end

Dist(Dist>R_m) = 0;

[Pair(:,2), Pair(:,1)] = find(Dist);

Pair = sortrows(sort(Pair,2));
Pair = Pair(1:2:end,:);

end

