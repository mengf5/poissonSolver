clear; close all;

% open file and interpret results

dir_open = './all_results/';
filename = 'results_20timesIterations_scanned.txt';

% open
fid = fopen([dir_open,filename],'r');

fgetl(fid);
[A,count] = fscanf(fid,'%d %d %d %d %d %d %e %e %e %e',[10,inf]);
data = A';

% Nx Ny Ngp nIter nRanks nThreads tInit tCalc tComm tTot
%  1  2   3     4      5        6     7     8     9   10

Nx = data(:,1);
Ny = data(:,2);
Ngp = data(:,3);
nIter = data(:,4);
Np = data(:,5); % nProcs
Nt = data(:,6); % nThreads
tInit = data(:,7);
tCalc = data(:,8);
tComm = data(:,9);
tTot = data(:,10);

% get averaged quantaties
tBarComm = tComm./(nIter);
tBarCalc = tCalc./(nIter.*Ngp);

for i = 1:length(data(:,1))
    if Nt(i) == 0
        Nt(i) = 1;
    end
end

figure(1)
plot((Nx.*Ny./(Nt.*Np) + Ngp.*Ny),tBarCalc,'*')
ylabel('tBarCalc')
figure(2)
plot(Np,tBarComm,'*')
ylabel('tBarComm')
% tBarCalc./(Nx.*Ny./(Nt.*Np) + Ngp.*Ny)
% tBarComm./(Np)


