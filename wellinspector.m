d = (-diff(wells)');

woi = 286;
% 20190701_water1_pump8_nocover
% 20190701_water2_pump8_nocover (exposure normalisazion)
% 20190621_lignin_pump8_nocover (autodetection doesn't work)
% 20190313_illite0_05
% 20190313_illite0_01
% 20190320_NaCl2_0
% 20190312_TapWater2 (good example of early onset detection in ROD)
% 20190313_illite0_1 (good example for missdetection of ROD)

figure
hold on
image(d)
scatter(Frz_Index,1:nWells,'xw')
scatter(Frz_Index_New,1:nWells,'+w')

ylabel('well nr')
% xlabel('T [°C]')
% xticklabels = (tdata(T_index(1):10:T_index(numel(T_index))));
% xticks = linspace(1, size(wells, 1), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xticks(1:numel(T_index));
xticklabels(pcrnm);
set(gca, 'TickLabelInterpreter', 'none')
xtickangle(90)
c = colorbar;
c.Label.String = 'diff(intensity)';

for i = 1:nTrays-1
yline(96*i,'--w','lineWidth',1)
end
ylim([1 nTrays*96])
xlim([min(Frz_Index)-1 max(Frz_Index)+1])

figure

hold on
xline(Frz_Index(woi),'--k')
xline(Frz_Index_New(woi),'-k')
plot((wells(:,woi)-mean(wells(:,woi)))/range((wells(:,woi)))+1)
plot(-diff(wells(:,woi))/range(-diff(wells(:,woi))))

legend('ROD','max(diff)')
