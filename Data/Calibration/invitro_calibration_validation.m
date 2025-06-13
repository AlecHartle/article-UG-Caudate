%% Upload Calibration Data for Plotting Figure 1h
clc, close all, clear all

%%%file path for invitro_calibration_table.mat
calibration_path = ''
fpath = calibration_path;
load(fpath);

%% parameters

marker_size = 10; % 12 is the default used in Seth's emostroop paper.
cap_size = 12;
line_width = 2;
xlimit  = [-200 2800];
ylimit = [-200 2800];
%% 

% Reference order, DA 5HT pH NE
% ------------------------------------- %
%              Graphing DA
% ------------------------------------- %
disp("Graphing DA")
DA_true = calibration_table.DA_table.DAMeanTrue;

DA_prediction_mean = calibration_table.DA_table.DAMeanPred;
DA_prediction_mean_SEM = calibration_table.DA_table.DASEM;

x5HT_prediction_mean = calibration_table.DA_table.x5HTResidualMeanPred;
x5HT_prediction_mean_SEM = calibration_table.DA_table.x5HTResidualSEM;

pH_prediction_mean = calibration_table.DA_table.pHResidualMeanPred;
pH_prediction_mean_SEM = calibration_table.DA_table.pHResidualSEM;

NE_prediction_mean = calibration_table.DA_table.NEResidualMeanPred;
NE_prediction_mean_SEM = calibration_table.DA_table.NEResidualSEM;


% fig_DA = figure;
figure;
hold on
e_pH = errorbar(DA_true, pH_prediction_mean, pH_prediction_mean_SEM,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',  [0.5000    0.5000    1.0000] ,'MarkerSize',marker_size, 'LineWidth', 1);
e_NE = errorbar(DA_true, NE_prediction_mean, NE_prediction_mean_SEM, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c','MarkerSize',marker_size, 'LineWidth', 1);
e_DA = errorbar(DA_true, DA_prediction_mean, DA_prediction_mean_SEM,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k','MarkerSize',marker_size, 'LineWidth', 1);
e_5HT = errorbar(DA_true, x5HT_prediction_mean, x5HT_prediction_mean_SEM, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm','MarkerSize',marker_size, 'LineWidth', 1);

DA_R = fitlm(DA_true,DA_prediction_mean);
DA_R_squared = DA_R.Rsquared.Ordinary;
disp(['DA R-squared: ', num2str(DA_R_squared)]);
plot(DA_true, DA_R.Fitted,'k', 'LineWidth', line_width);


e_DA.CapSize = cap_size;
e_DA.Color = 'k';
e_5HT.CapSize = cap_size;
e_5HT.Color = 'k';
e_pH.CapSize = cap_size;
e_pH.Color = 'k';
e_NE.CapSize = cap_size;
e_NE.Color = 'k';

ylim(ylimit);
xlim(xlimit);
xticks([0 500 1000 1500 2000 2500]);
yticks([0 500 1000 1500 2000 2500]);
xlabel(sprintf('True [ DA or 5HT or NE or pH ] (nM)'));
ylabel(sprintf('Predicted [ DA ] (nM)'));

ax = gca;
ax.FontSize = 18;
axis square
box on
boxLineWidth = 3;
set(gca, 'LineWidth', boxLineWidth);

hold off

%% 

% ------------------------------------- %
%             Graphing NE
% ------------------------------------- %
disp("Graphing NE")
NE_true = calibration_table.NE_table.NEMeanTrue;

DA_prediction_mean = calibration_table.NE_table.DAResidualMeanPred;
DA_prediction_mean_SEM = calibration_table.NE_table.DAResidualSEM;

x5HT_prediction_mean = calibration_table.NE_table.x5HTResidualMeanPred;
x5HT_prediction_mean_SEM = calibration_table.NE_table.x5HTResidualSEM;

pH_prediction_mean = calibration_table.NE_table.pHResidualMeanPred;
pH_prediction_mean_SEM = calibration_table.NE_table.pHResidualSEM;

NE_prediction_mean = calibration_table.NE_table.NEMeanPred;
NE_prediction_mean_SEM = calibration_table.NE_table.NESEM;

fig_NE = figure;
hold on
e_pH = errorbar(NE_true, pH_prediction_mean, pH_prediction_mean_SEM,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',  [0.5000    0.5000    1.0000] ,'MarkerSize',marker_size, 'LineWidth', 1);
e_NE = errorbar(NE_true, NE_prediction_mean, NE_prediction_mean_SEM, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c','MarkerSize',marker_size, 'LineWidth', 1);
e_DA = errorbar(NE_true, DA_prediction_mean, DA_prediction_mean_SEM,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k','MarkerSize',marker_size, 'LineWidth', 1);
e_5HT = errorbar(NE_true, x5HT_prediction_mean, x5HT_prediction_mean_SEM, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm','MarkerSize',marker_size, 'LineWidth', 1);

NE_R = fitlm(NE_true, NE_prediction_mean);
NE_R_squared = NE_R.Rsquared.Ordinary;
disp(['NE R-squared: ', num2str(NE_R_squared)]);
plot(NE_true, NE_R.Fitted,'k', 'LineWidth', line_width);
hold off

e_DA.CapSize = cap_size;
e_DA.Color = 'k';
e_5HT.CapSize = cap_size;
e_5HT.Color = 'k';
e_pH.CapSize = cap_size;
e_pH.Color = 'k';
e_NE.CapSize = cap_size;
e_NE.Color = 'k';

ylim(ylimit);
xlim(xlimit);
xticks([0 500 1000 1500 2000 2500]);
yticks([0 500 1000 1500 2000 2500]);
xlabel(sprintf('True [ DA or 5HT or NE or pH ] (nM)'));
ylabel(sprintf('Predicted [ NE ] (nM)'));

ax = gca;
ax.FontSize = 18;
axis square
box on
boxLineWidth = 3;
set(gca, 'LineWidth', boxLineWidth);


%% 
% ------------------------------------- %
%             Graphing 5HT
% ------------------------------------- %
disp("Graphing 5HT")
x5HT_true = calibration_table.x5HT_table.x5HTMeanTrue;

DA_prediction_mean = calibration_table.x5HT_table.DAResidualMeanPred;
DA_prediction_mean_SEM = calibration_table.x5HT_table.DAResidualSEM;

x5HT_prediction_mean = calibration_table.x5HT_table.x5HTMeanPred;
x5HT_prediction_mean_SEM = calibration_table.x5HT_table.x5HTSEM;

pH_prediction_mean = calibration_table.x5HT_table.pHResidualMeanPred;
pH_prediction_mean_SEM = calibration_table.x5HT_table.pHResidualSEM;

NE_prediction_mean = calibration_table.x5HT_table.NEResidualMeanPred;
NE_prediction_mean_SEM = calibration_table.x5HT_table.NEResidualSEM;


fig_5HT = figure;
hold on
e_pH = errorbar(x5HT_true, pH_prediction_mean, pH_prediction_mean_SEM,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',  [0.5000    0.5000    1.0000] ,'MarkerSize',marker_size, 'LineWidth', 1);
e_NE = errorbar(x5HT_true, NE_prediction_mean, NE_prediction_mean_SEM, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c','MarkerSize',marker_size, 'LineWidth', 1);
e_DA = errorbar(x5HT_true, DA_prediction_mean, DA_prediction_mean_SEM,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k','MarkerSize',marker_size, 'LineWidth', 1);
e_5HT = errorbar(x5HT_true, x5HT_prediction_mean, x5HT_prediction_mean_SEM, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm','MarkerSize',marker_size, 'LineWidth', 1);

x5HT_R = fitlm(x5HT_true,x5HT_prediction_mean);
x5HT_R_squared = x5HT_R.Rsquared.Ordinary;
disp(['x5HT R-squared: ', num2str(x5HT_R_squared)]);
plot(x5HT_true, x5HT_R.Fitted,'k', 'LineWidth', line_width);
hold off

e_DA.CapSize = cap_size;
e_DA.Color = 'k';
e_5HT.CapSize = cap_size;
e_5HT.Color = 'k';
e_pH.CapSize = cap_size;
e_pH.Color = 'k';
e_NE.CapSize = cap_size;
e_NE.Color = 'k';

ylim(ylimit);
xlim(xlimit);
xticks([0 500 1000 1500 2000 2500]);
yticks([0 500 1000 1500 2000 2500]);
xlabel(sprintf('True [ DA or 5HT or NE or pH ] (nM)'));
ylabel(sprintf('Predicted [ 5HT ] (nM)'));

ax = gca;
ax.FontSize = 18;
axis square
box on
boxLineWidth = 3;
set(gca, 'LineWidth', boxLineWidth);


