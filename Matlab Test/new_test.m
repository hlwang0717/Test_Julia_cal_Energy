
info = load('jamming_config.txt');
Config = get_RandomConfig(0.84, 1024, 2);
Config.PositX = info(:,1);
Config.PositY = info(:,2);
Config.Rad = info(:,3);
Config.Type = info(:,4);
Config.Phi = pi*sum(Config.Rad.^2);

tic
Pair = find_Pair(Config);
nlist = NeighbourList(Config, 5); NPair = length(nlist.i);
x = [Config.PositX; Config.PositY];

for i = 1:1000
    [E,~] = cal_Energy(x, Config, Pair, NPair);
end
time1 = toc

% tic
% for i = 1:500
%     nlist = NeighborList(Config, 2.5);
% end
% time2 = toc

% tic
% Config1 = min_Energy(Config);
% tim2 = toc


