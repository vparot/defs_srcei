% Analisis de defunciones inscritas en el Registro Civil.
% Genera graficos separados para todo Chile y para cada region. Cada
% grafico muestra una estimacion de las defunciones por semana basada en
% proyeccion log-lineal entre los años 2010-2019. La desviacion estandar
% del error de estimacion por semana se muestra en el grafico sobre y bajo
% la curva estimada. Las defunciones inscritas en el año 2020 se muestran
% directamente en el grafico. 
%
% 2020 Vicente Parot
% Wellman Center for Photomedicine
% Harvard Medical School
%
% Instituto de Ingeniería Biológica y Médica
% Pontificia Universidad Católica de Chile
%

%% calculate stats
load defunciones
defs.DayOfYear = day(datetime(defs{:,1},defs{:,2},defs{:,3}),'dayofyear');
defs.WeekOfYear = week(datetime(defs{:,1},defs{:,2},defs{:,3}));

uAnyos = unique(defs.A_O);
uComuna = unique(defs.COMUNA);
uRegion = unique(defs.REGION);

%% show data per week, all regions
% collect data
allwks = [];
for it = 1:numel(uAnyos)
    mAnyo = defs.A_O == uAnyos(it);
    mSel = mAnyo;
    totPerWkThisYr = varfun(@sum,defs(mSel,:),'InputVariables','TOTAL','GroupingVariables','WeekOfYear');
    totals = totPerWkThisYr.sum_TOTAL;
    totals(end:53) = nan;
    allwks = [allwks totals];
end
% collect first and last week
allwks(1) = allwks(2);
allwks(1,2:end) = allwks(1,2:end) + nansum(allwks(52:end,1:end-1));
allwks(52:end,:) = [];
% median filtering to minimize weekend and holiday artifacts
allwks = medfilt2(allwks,[5 1]);
% extract trend from first principal component
[u, s, v] = svd(allwks(1:end,1:end-1));
% fit log-linear increase in rates
projmat = ((1:size(v,1)+1)'*[0 1]+[1 0]);
coeffs = projmat(1:end-1,:)\log(v(:,1));
sel = 1:1;
allwksest = abs(u(:,sel)*s(sel,sel)*exp(projmat*coeffs)');

% make figure
figure
xax = (1:size(allwks,1))';
ydat = allwksest(:,end) + [-1 1].*std(allwks(:,1:end-1)-allwksest(:,1:end-1),[],2);%/sqrt(size(allwks,2)-1);
xax = xax(:);
ydat = ydat(:);
cli = colormap(lines);
plot(xax,[allwksest(:,end)],'color',cli(1,:),'linewidth',2)
hold on
plot(xax,[allwks(:,end)],'color',cli(2,:),'linewidth',2)
patch([xax; xax(end:-1:1)],[ydat(1:end/2); ydat(end:-1:end/2+1)],cli(1,:),'edgecolor','none','facealpha',.2)
ylim(ylim.*[0 1])
xlim([1 51])
xlabel 'Semana'
ylabel 'Defunciones'
legend({'Proyección ajustada desde 2010',['Defunciones 2020 (hasta ' fecha ')']},'location','se')
title({'Defunciones por cualquier causa registradas en Chile','fuente: estadísticas Registro Civil'})
saveas(gcf,'out_00.png')    

%% show data per week, each region
for itReg = 1:numel(uRegion)-1
    if strcmp(uRegion{itReg},{'XVI Región del Ñuble'})
        continue
    end
    allwks = [];
    for it = 1:numel(uAnyos)
        mAnyo = defs.A_O == uAnyos(it);
        mRegion = strcmp(defs.REGION,uRegion{itReg});
        if strcmp(uRegion{itReg},{'VIII Región de Concepción'})
            mRegion = mRegion | strcmp(defs.REGION,'XVI Región del Ñuble');
        end
        mSel = mAnyo & mRegion;
        totPerWkThisYr = varfun(@sum,defs(mSel,:),'InputVariables','TOTAL','GroupingVariables','WeekOfYear');
        totals = totPerWkThisYr.sum_TOTAL;
        totals(end:53) = nan;
        allwks = [allwks totals];
    end
    % collect first and last week
    allwks(1) = allwks(2);
    allwks(1,2:end) = allwks(1,2:end) + nansum(allwks(52:end,1:end-1));
    allwks(52:end,:) = [];
    % median filtering to minimize weekend and holiday artifacts
    allwks = medfilt2(allwks,[5 1]);
    % extract trend from first principal component
    [u, s, v] = svd(allwks(1:end,1:end-1));
    % fit log-linear increase in rates
    projmat = ((1:size(v,1)+1)'*[0 1]+[1 0]);
    coeffs = projmat(1:end-1,:)\log(v(:,1));
    sel = 1:1;
    allwksest = abs(u(:,sel)*s(sel,sel)*exp(projmat*coeffs)');

    % make figure
    figure
    xax = (1:size(allwks,1))';
    ydat = allwksest(:,end) + [-1 1].*std(allwks(:,1:end-1)-allwksest(:,1:end-1),[],2);%/sqrt(size(allwks,2)-1);
    xax = xax(:);
    ydat = ydat(:);
    cli = colormap(lines);
    plot(xax,[allwksest(:,end)],'color',cli(1,:),'linewidth',2)
    hold on
    plot(xax,[allwks(:,end)],'color',cli(2,:),'linewidth',2)
    patch([xax; xax(end:-1:1)],[ydat(1:end/2); ydat(end:-1:end/2+1)],cli(1,:),'edgecolor','none','facealpha',.2)
    ylim(ylim.*[0 1])
    xlim([1 51])
    xlabel 'Semana'
    ylabel 'Defunciones'
    legend({'Proyección ajustada desde 2010',['Defunciones 2020 (hasta ' fecha ')']},'location','se')
    title({'Defunciones por cualquier causa',uRegion{itReg}})
    saveas(gcf,sprintf('out_%02d.png',itReg))
end

%% datos por semana, por comuna, en regiones seleccionadas
% for itReg = [6 7 13]
%     if strcmp(uRegion{itReg},{'XVI Región del Ñuble'})
%         continue
%     end
%     mRegion = strcmp(defs.REGION,uRegion{itReg});
%     uComunaReg = unique(defs(mRegion,:).COMUNA);
%     for itCom = 1:numel(uComunaReg)
%         mComuna = strcmp(defs.COMUNA,uComunaReg{itCom});
%         allwks = [];
%         for it = 1:numel(uAnyos)
%             mAnyo = defs.A_O == uAnyos(it);
%             if strcmp(uRegion{itReg},{'VIII Región de Concepción'})
%                 mRegion = mRegion | strcmp(defs.REGION,'XVI Región del Ñuble');
%             end
%             mSel = mAnyo & mRegion & mComuna;
%             totPerWkThisYr = varfun(@sum,defs(mSel,:),'InputVariables','TOTAL','GroupingVariables','WeekOfYear');
%             totals = zeros(53,1);
%             totals(totPerWkThisYr.WeekOfYear) = totals(totPerWkThisYr.WeekOfYear) + totPerWkThisYr.sum_TOTAL;
%             if uAnyos(it) == 2020
%                 thisWk = find(totPerWkThisYr.WeekOfYear==max(totPerWkThisYr.WeekOfYear));
%                 if isempty(thisWk)
%                     % algunas comunas como juan fernandez y primavera no tienen muertes todos los años.
%                     disp(defs(mRegion & mComuna,:))
%                     continue
%                 end
%                 totals(thisWk) = totals(thisWk)*max(totPerWkThisYr.GroupCount)/totPerWkThisYr(thisWk,:).GroupCount;
%                 totals(thisWk+1:end) = nan;
%             end
%             allwks = [allwks totals]; 
%         end
%         if any(sum(~~allwks(:,1:end-1))/53 < .9)
%             continue
%         end
%         avgHistPerWeek = mean(sum(allwks(:,1:end-1)))/52;
%         if avgHistPerWeek < 1
%             continue
%         end
%         sum(allwks);
%         % collect first and last week
% %         allwks(1) = allwks(2); 
%         allwks(1,2:end) = allwks(1,2:end) + allwks(53:end,1:end-1);
%         allwks(53:end,:) = [];
%         % median filtering to minimize weekend and holiday artifacts
%         allwks(1:thisWk,end) = reducepeaks(allwks(1:thisWk,end));
%         allwks(1:end,1:end-1) = reducepeaks(allwks(1:end,1:end-1));
%         % extract trend from first principal component
%         [u, s, v] = svd(allwks(1:end,1:end-1));
%         % fit log-linear increase in rates
%         projmat = ((1:size(v,1)+1)'*[0 1]+[1 0]);
%         coeffs = projmat(1:end-1,:)\log(v(:,1));
%         sel = 1:1;
%         allwksest = abs(u(:,sel)*s(sel,sel)*exp(projmat*coeffs)');
%  
%         % make figure
%         figure windowstyle docked
%         xax = (1:size(allwks,1))';
%         ydat = allwksest(:,end) + [-1 1].*std(allwks(:,1:end-1)-allwksest(:,1:end-1),[],2);%/sqrt(size(allwks,2)-1);
%         xax = xax(:);
%         ydat = ydat(:);
%         cli = colormap(lines);
%         plot(xax,[allwksest(:,end)],'color',cli(1,:),'linewidth',2)
%         hold on
%         plot(xax,[allwks(:,end)],'color',cli(2,:),'linewidth',2)
%         patch([xax; xax(end:-1:1)],[ydat(1:end/2); ydat(end:-1:end/2+1)],cli(1,:),'edgecolor','none','facealpha',.2)
%         ylim(ylim.*[0 1])
%         xlim([1 51])
%         xlabel 'Semana'
%         ylabel 'Defunciones'
%         legend({'Proyección ajustada desde 2010',['Defunciones 2020 (hasta ' fecha ')']},'location','se')
%         title({'Defunciones por cualquier causa',[uRegion{itReg} ', comuna de ' capitalize(uComunaReg{itCom})]})
%     %     saveas(gcf,sprintf('out_%02d.png',itReg))
%     end % for Comuna
% end % for Region
