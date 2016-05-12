
clear; close all;

% open file and interpret results

fig = 1;

dir_open = './all_results/';

nodes = [8,16,32];
for i = 1:length(nodes)
    node = nodes(i);
    filename = ['newest_results_',num2str(node),'_scanned.txt'];

    % open
    fid = fopen([dir_open,filename],'r');

    fgetl(fid);
    [A,count] = fscanf(fid,'%d %d %d %d %d %d %e %e %e %e',[10,inf]);
    data = A';

    % Nx Ny Ngp nIter nRanks nThreads tInit tCalc tComm tTot
    %  1  2   3     4      5        6     7     8     9   10

    % put data into vecors
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

    % data for each ngp config
    data1 = [];
    data2 = [];
    data4 = [];
    data8 = [];
    
    % collect data
    for j = 1:length(Ngp)
        if Ngp(j) == 1
            data1 = [data1;Np(j),tTot(j),tComm(j),tCalc(j)];
        elseif Ngp(j) == 2
            data2 = [data2;Np(j),tTot(j),tComm(j),tCalc(j)];
        elseif Ngp(j) == 4
            data4 = [data4;Np(j),tTot(j),tComm(j),tCalc(j)];
        elseif Ngp(j) == 8
            data8 = [data8;Np(j),tTot(j),tComm(j),tCalc(j)];
        else
            error('unexpected value for Ngp')
        end
    end
    
    % sort data
    [~,I] = sort(data1(:,1),1,'ascend');
    data1 = data1(I,:);
    [~,I] = sort(data2(:,1),1,'ascend');
    data2 = data2(I,:);
    [~,I] = sort(data4(:,1),1,'ascend');
    data4 = data4(I,:);
    [~,I] = sort(data8(:,1),1,'ascend');
    data8 = data8(I,:);
    

    figure(fig); hold on
    fig = fig+1;
    plot(data1(:,1),data1(:,2),'r*--')
    plot(data2(:,1),data2(:,2),'g*--')
    plot(data4(:,1),data4(:,2),'b*--')
    plot(data8(:,1),data8(:,2),'m*--')
    legend({'N_{gp} = 1','N_{gp} = 2','N_{gp} = 4','N_{gp} = 8'})
    xlabel('Number of Processors')
    ylabel('Execution Time (s)')
    title([num2str(node),' nodes'])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'Xtick',[data1(:,1)])
    grid on
    
    
    fclose(fid);
end
