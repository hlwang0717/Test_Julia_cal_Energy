function config = get_RandomConfig( Phi, Num, Alpha )

Length = 1.0; 
Height = 1.0;

BC = 0;

StrainType = 1;
Strain = 0;

E0 = 1.0;

Particle_Num = Num;

X = Length*rand(Particle_Num,1);
Y = Height*rand(Particle_Num,1);

r_ratio = 1.4;
r = sqrt(2*Phi*Length*Height/pi/Particle_Num/(1+r_ratio^2));

Id = (1:Particle_Num)';
Type = randerr(1, Particle_Num, [0, Particle_Num/2; 0, 1])';

Rad = ((r_ratio-1)*Type + 1)*r;

config.BC = BC;
config.Length = Length;
config.Height = Height;
config.Num = Particle_Num;
config.E0 = E0;
config.Alpha = Alpha;
config.Id = Id;
config.Type = Type;
config.PositX = X;
config.PositY = Y;
config.Rad = Rad;
config.Phi = Phi;

config.StrainType = StrainType;
config.Strain = Strain;

config.ref = config;

% Pair = find_Pair(config); NPair = length(Pair(:,1));
% config.Energy = cal_Energy([X;Y], config, Pair, NPair);
% config.Pressure = cal_Pressure([X;Y], config, Pair);

% % while true
% %     deform.type = 'grow';
% %     deform.value = 1;
% %     config = deform_Config(config, deform);
% %     if config.Energy < 1.0001*eps
% %         break
% %     end
% % end

end

