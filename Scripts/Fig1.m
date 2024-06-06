% nanmedian([TwoTapmsfinaldata.Age, TwoTapmsfinaldata.sesscore, TwoTapmsfinaldata.Mindfulness, TwoTapmsfinaldata.Depression])
% mad([TwoTapmsfinaldata.Age, TwoTapmsfinaldata.sesscore, TwoTapmsfinaldata.Mindfulness, TwoTapmsfinaldata.Depression])
% length(find(TwoTapmsfinaldata.gender==1))
% length(find(TwoTapmsfinaldata.gender==2))
% length(find(TwoTapmsfinaldata.gender==3))
% length(find(TwoTapmsfinaldata.gender==1))./length(find(TwoTapmsfinaldata.gender))
% length(find(TwoTapmsfinaldata.gender==2))./length(find(TwoTapmsfinaldata.gender))
% length(find(TwoTapmsfinaldata.gender==3))./length(find(TwoTapmsfinaldata.gender))
% length(find(TwoTapmsfinaldata.ethnicity==1))
% length(find(TwoTapmsfinaldata.ethnicity==2))
% length(find(TwoTapmsfinaldata.ethnicity==3))
% length(find(TwoTapmsfinaldata.ethnicity==4))
% length(find(TwoTapmsfinaldata.ethnicity==5))
% length(find(TwoTapmsfinaldata.ethnicity==6))
% 
% length(find(TwoTapmsfinaldata.ethnicity==1))./length(find(TwoTapmsfinaldata.ethnicity))
% length(find(TwoTapmsfinaldata.ethnicity==2))./length(find(TwoTapmsfinaldata.ethnicity))
% length(find(TwoTapmsfinaldata.ethnicity==3))./length(find(TwoTapmsfinaldata.ethnicity))
% length(find(TwoTapmsfinaldata.ethnicity==4))./length(find(TwoTapmsfinaldata.ethnicity))
% length(find(TwoTapmsfinaldata.ethnicity==5))./length(find(TwoTapmsfinaldata.ethnicity))
% length(find(TwoTapmsfinaldata.ethnicity==6))./length(find(TwoTapmsfinaldata.ethnicity))

data = readtable('A:\TwoTap\Manuscript\TwoTap_ms_finaldata_May29.xlsx');
data = renamevars(data,["age"],["Age"]);

data.gw_efficiency_go = str2double(data.gw_efficiency_go);
outcon = find(isoutlier(data.Consistency));
data.Consistency(outcon) = nan;
%length(find(isnan(data.Consistency)));

zCon = (data.Consistency - nanmean(data.Consistency))./nanstd(data.Consistency);
data.zCon = zCon;

zAge = (data.Age - nanmean(data.Age))./nanstd(data.Age);
data.zAge = (zAge);

zgw = (data.gw_efficiency_go - nanmean(data.gw_efficiency_go))./nanstd(data.gw_efficiency_go);
data.zgw = (zgw);

zmf = (data.mf_efficiency - nanmean(data.mf_efficiency))./nanstd(data.mf_efficiency);
data.zmf = (zmf);

zls = (data.ls_span_efficiency - nanmean(data.ls_span_efficiency))./nanstd(data.ls_span_efficiency);
data.zls = (zls);

zfo = (data.fo_efficiency - nanmean(data.fo_efficiency))./nanstd(data.fo_efficiency);
data.zfo = (zfo);

zglobal = (data.global_efficiency - nanmean(data.global_efficiency))./nanstd(data.global_efficiency);
data.zglobal = (zglobal);

subsetData = data(data.Age>60,:);

%mdlspec = 'Consistency ~ Age + gender + ethnicity + sesscore + Mindfulness + Depression';
mdlspec = 'zCon ~ zAge + gender + ethnicity + sesscore + Mindfulness + Depression';
% mdlspec = 'zCon ~ AgeBin + gender + ethnicity + sesscore + Mindfulness + Depression';

model = fitlm(data,mdlspec,'RobustOpts','on'); %'CategoricalVars', [3 4]
%zAge             -0.23802    0.088037     -2.7036    0.007517
model = fitlm(subsetData,mdlspec,'RobustOpts','on') %'CategoricalVars', [3 4]
%zAge            -0.077083     0.31938    -0.24135    0.80964
% hist(data.Consistency)

%%% Inhibitory Control 
mdlspec = 'zgw ~ zCon + Age';
model = fitlm(data,mdlspec,'RobustOpts','on')
%zCon             0.18137     0.062764     2.8897     0.0044131

%%% Interference Processing
mdlspec = 'zmf ~ zCon + Age';
model = fitlm(data,mdlspec,'RobustOpts','on')
%zCon           0.073357     0.058928     1.2449       0.21426

%%% Working Memory
mdlspec = 'zls ~ zCon + Age';
model = fitlm(data,mdlspec,'RobustOpts','on')
%zCon            0.089428     0.059874     1.4936       0.13643

%%% Emotion Bias
mdlspec = 'zfo ~ zCon + Age';
model = fitlm(data,mdlspec,'RobustOpts','on')
%zCon            0.055175     0.054881     1.0054        0.3156

%%% Global
mdlspec = 'zglobal ~ zCon + Age';
model = fitlm(data,mdlspec,'RobustOpts','on')
%zCon              0.1331     0.053043     2.5094      0.012666

% plot(data.Consistency, data.gw_efficiency_go, '.')
% plot(data.Consistency, data.mf_efficiency, '.')
% plot(data.Consistency, data.ls_span_efficiency, '.')					
% plot(data.Consistency, data.fo_efficiency, '.')	
% plot(data.Consistency, data.global_efficiency, '.')



%% 1D
close all
figure(4)
X=data.Consistency;
Y=data.Age;

badSubs=isnan(X) | isnan(Y);
coefficients = polyfit(X(~badSubs), Y(~badSubs), 1);
xFit = linspace(min(X), max(X), 1000);
yFit = polyval(coefficients , xFit);
hold on
plot(X, Y,'.','MarkerSize',8)
plot(xFit, yFit, 'r-', 'LineWidth', 1.5); % Plot fitted line.
ylabel('Age (years)','fontsize',12)
xlabel('Interoceptive Consistency','fontsize',12)
f=gcf;
f.Position = [1000 918 560 420];
exportgraphics(gcf, "A:\TwoTap\Manuscript\Figure\Fig1D.tif" )
%% 1E
close all

figure(5)
X=data.Consistency;
Y=data.global_efficiency;

badSubs=isnan(X) | isnan(Y);
coefficients = polyfit(X(~badSubs), Y(~badSubs), 1);
xFit = linspace(min(X), max(X), 1000);
yFit = polyval(coefficients , xFit);
hold on
plot(X, Y,'.','MarkerSize',8)
plot(xFit, yFit, 'r-', 'LineWidth', 1.5); % Plot fitted line.
ylabel('Global Cognitive Efficiency','fontsize',12)
xlabel('Interoceptive Consistency','fontsize',12)
f=gcf;
f.Position = [1000 918 560 420];
exportgraphics(gcf, "A:\TwoTap\Manuscript\Figure\Fig1E.tif" )

%% Fig 1E-2

figure(6)
subplot(2,2,1)

X=data.Consistency;
Y=data.gw_efficiency_go;
badSubs=isnan(X) | isnan(Y);
coefficients = polyfit(X(~badSubs), Y(~badSubs), 1);
xFit = linspace(min(X), max(X), 1000);
yFit = polyval(coefficients , xFit);
hold on
plot(X, Y,'.','MarkerSize',8)
plot(xFit, yFit, 'r-', 'LineWidth', 1.5); % Plot fitted line.
ylabel('Inhibitory Control Eff.','fontsize',12)
xlabel('Interoceptive Consistency','fontsize',12)

subplot(2,2,2)
X=data.Consistency;
Y=data.ls_span_efficiency;
badSubs=isnan(X) | isnan(Y);
coefficients = polyfit(X(~badSubs), Y(~badSubs), 1);
xFit = linspace(min(X), max(X), 1000);
yFit = polyval(coefficients , xFit);
hold on
plot(X, Y,'.','MarkerSize',8)
plot(xFit, yFit, 'r-', 'LineWidth', 1.5); % Plot fitted line.
ylabel('Working Memory Eff.','fontsize',12)
xlabel('Interoceptive Consistency','fontsize',12)

subplot(2,2,3)
X=data.Consistency;
Y=data.fo_efficiency;
badSubs=isnan(X) | isnan(Y);
coefficients = polyfit(X(~badSubs), Y(~badSubs), 1);
xFit = linspace(min(X), max(X), 1000);
yFit = polyval(coefficients , xFit);
hold on
plot(X, Y,'.','MarkerSize',8)
plot(xFit, yFit, 'r-', 'LineWidth', 1.5); % Plot fitted line.
ylabel('Emotion Bias Eff.','fontsize',12)
xlabel('Interoceptive Consistency','fontsize',12)

subplot(2,2,4)
X=data.Consistency;
Y=data.mf_efficiency;
badSubs=isnan(X) | isnan(Y);
coefficients = polyfit(X(~badSubs), Y(~badSubs), 1);
xFit = linspace(min(X), max(X), 1000);
yFit = polyval(coefficients , xFit);
hold on
plot(X, Y,'.','MarkerSize',8)
plot(xFit, yFit, 'r-', 'LineWidth', 1.5); % Plot fitted line.
ylabel('Interference Proc. Eff.','fontsize',12)
xlabel('Interoceptive Consistency','fontsize',12)
exportgraphics(gcf, "A:\TwoTap\Manuscript\Figure\Fig1E2.tif" )



