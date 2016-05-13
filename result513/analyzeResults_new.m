
clear; close all;

% open file and interpret results

fig = 1;

dir_open = './';

nodes = [8,16,32];
for i = 1:length(nodes)
    node = nodes(i);
    filename = ['newest_results_5_13_',num2str(node),'_advanced.txt'];

    % open
    fid = fopen([dir_open,filename],'r');

    fgetl(fid);
    [A,count] = fscanf(fid,'%d %d %d %d %d %d %e %e %e %e %e %e %e',[13,inf]);
    data = A';

    % Nx    Ny   Ngp  nIter nRanks  nThreads      tInit      tCalc      
    %  1     2     3      4      5         6          7          8
    % tComm     tbatch     tTotal   tC/tcal    (tC+tb)/tCal
    %     9         10         11        12              13 
    
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
    tBatch = data(:,10);
    tTot = data(:,11);
    ratio1 = data(:,12);
    ratio2 = data(:,13);

    % data for each ngp config
    data1 = [];
    data16 = [];
    data4 = [];
    data8 = [];
    
    % collect data
    for j = 1:length(Ngp)
        if Ngp(j) == 1
            data1 = [data1;Np(j),tTot(j),ratio1(j),ratio2(j)];
        elseif Ngp(j) == 16
            data16 = [data16;Np(j),tTot(j),ratio1(j),ratio2(j)];
        elseif Ngp(j) == 4
            data4 = [data4;Np(j),tTot(j),ratio1(j),ratio2(j)];
        elseif Ngp(j) == 8
            data8 = [data8;Np(j),tTot(j),ratio1(j),ratio2(j)];
        else
            Ngp(j)
            error('unexpected value for Ngp')
        end
    end
    
    % sort data
    [~,I] = sort(data1(:,1),1,'ascend');
    data1 = data1(I,:);
    [~,I] = sort(data4(:,1),1,'ascend');
    data4 = data4(I,:);
    [~,I] = sort(data8(:,1),1,'ascend');
    data8 = data8(I,:);
    [~,I] = sort(data16(:,1),1,'ascend');
    data16 = data16(I,:);
    
    % plot execution time/iter
    figure(fig); hold on
    fig = fig+1;
    plot(data1(:,1),data1(:,2),'r*--')
    plot(data4(:,1),data4(:,2),'g*--')
    plot(data8(:,1),data8(:,2),'b*--')
    plot(data16(:,1),data16(:,2),'m*--')
    legend({'N_{gp} = 1','N_{gp} = 4','N_{gp} = 8','N_{gp} = 16'})
    xlabel('Number of Processors')
    ylabel('t_{tot} / N_{tot}   (sec / iter)')
    title([num2str(node),' nodes'])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'Xtick',[data1(:,1)])
    set(gca,'Color','White')
    xlim([data1(1,1)*.5,data1(end,1)*4])
    grid on
    set(gca,'FontSize',12)
    xl = get(gca,'xlabel');
    set(xl,'FontSize',12)
    yl = get(gca,'ylabel');
    set(yl,'FontSize',12)
    tl = get(gca,'title');
    set(tl,'FontSize',12)
    grid on
    set(gcf,'Color','w')

    % plot tC/tCalc
    figure(fig); hold on
    fig = fig+1;
    plot(data1(:,1),data1(:,3),'r*--')
    plot(data4(:,1),data4(:,3),'g*--')
    plot(data8(:,1),data8(:,3),'b*--')
    plot(data16(:,1),data16(:,3),'m*--')
    legend({'N_{gp} = 1','N_{gp} = 4','N_{gp} = 8','N_{gp} = 16'})
    xlabel('Number of Processors')
    ylabel('t_{comm} / t_{calc}')
    title([num2str(node),' nodes'])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'Xtick',[data1(:,1)])
    set(gca,'Color','White')
    xlim([data1(1,1)*.5,data1(end,1)*4])
    grid on
    set(gca,'FontSize',12)
    xl = get(gca,'xlabel');
    set(xl,'FontSize',12)
    yl = get(gca,'ylabel');
    set(yl,'FontSize',12)
    tl = get(gca,'title');
    set(tl,'FontSize',12)
    grid on
    set(gcf,'Color','w')
    
    % plot (tC+tB)/tCalc
    figure(fig); hold on
    fig = fig+1;
    plot(data1(:,1),data1(:,4),'r*--')
    plot(data4(:,1),data4(:,4),'g*--')
    plot(data8(:,1),data8(:,4),'b*--')
    plot(data16(:,1),data16(:,4),'m*--')
    legend({'N_{gp} = 1','N_{gp} = 4','N_{gp} = 8','N_{gp} = 16'})
    xlabel('Number of Processors')
    ylabel(['(t_{comm}+t_{batch}) / t_{calc} '])
    title([num2str(node),' nodes'])
    xlim([data1(1,1)*.5,data1(end,1)*4])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'Xtick',[data1(:,1)])
    set(gca,'Color','White')
    set(gca,'FontSize',12)
    xl = get(gca,'xlabel');
    set(xl,'FontSize',12);
    yl = get(gca,'ylabel');
    set(yl,'FontSize',12);
    tl = get(gca,'title');
    set(tl,'FontSize',12);
    grid on
    set(gcf,'Color','w')
    
    fclose(fid);
end
