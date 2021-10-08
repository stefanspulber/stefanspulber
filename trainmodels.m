%% train and test multiple regression model 
clc, tic, disp('Working...')
varidx=[4 7 9 11:21];
for idx=1:6
%     mdl=stepwiselm(train, 'ResponseVar', 'madrs1', 'PredictorVars', varidx, 'Criterion', 'sse', 'PEnter', 0.09,'NSteps', idx, 'Lower', 'constant', 'Upper', 'linear', 'Verbose', 0);
    mdl=stepwiselm(train, 'ResponseVar', 'madrs1', 'PredictorVars', varidx, 'Criterion', 'aic', 'NSteps', idx, 'Lower', 'constant', 'Upper', 'linear', 'Verbose', 0);
    
    % collect info on trained model
    mdlh(idx).varsin=mdl.PredictorNames; 
    mdlh(idx).weights=mdl.Coefficients;
    mdlh(idx).train=[train.madrs1 mdl.Fitted]; 
    mdlh(idx).performance=[mdl.Rsquared.Adjusted mdl.RMSE];
    r0=corrcoef(table2array(train(:, mdl.PredictorNames)));
    mdlh(idx).vif=diag(inv(r0)); % VIF in training dataset
    
    % test performance on test dataset
    [ypred, yci]=predict(mdl, test(:, varidx));
    mdlh(idx).test=[test.madrs1 ypred yci];
    r=corrcoef(test.madrs1, ypred);
    r0=corrcoef(table2array(test(:, mdl.PredictorNames)));
    mdlh(idx).testp=[r(2)^2 sqrt(mean((test.madrs1-ypred).^2))];
    mdlh(idx).tvif=diag(inv(r0)); % VIF in test dataset
    mdlh(idx).precision=[sum(abs(diff(mdlh(idx).test(:, 1:2), 1, 2))<1); sum(abs(diff(mdlh(idx).test(:, 1:2), 1, 2))<2);...
        sum(abs(diff(mdlh(idx).test(:, 1:2), 1, 2))<3); sum(abs(diff(mdlh(idx).test(:, 1:2), 1, 2))<4);...
        sum(abs(diff(mdlh(idx).test(:, 1:2), 1, 2))<5)]/12;
end
toc

a=reshape([mdlh.performance], 2, [])';
b=reshape([mdlh.testp], 2, [])';
figure, plot(a(:, 2), a(:, 1), 'bo', b(:, 2), b(:, 1), 'ro'); grid on;
set(gca, 'ylim', [0 1]); title('train -> test');


%% plot correlation matrices
% varidx=[5 13 14 19 20 8 10 12 7 4]; % select only predictors included in models - sorted logically
varidx=[5 6 4 9:23 34:35 25:26 28:30]; % all predictors in training
% varidx=[5 15 16 21 22 10 12 14 9 4 29]; % select predictors for madrs12
% models
lbls=trainf.Properties.VariableNames(varidx); % get predictor names

lim=1; ramp=100;
corrcmap=[zeros(1, lim) linspace(0, 1, ramp) ones(1, 256-lim-ramp);...
    zeros(1, lim) linspace(0, 1, ramp) ones(1, 256-2*lim-2*ramp) linspace(1, 0, ramp) zeros(1, lim);...
    ones(1, 256-lim-ramp) linspace(1, 0, ramp) zeros(1, lim)];

figure, 
% trainmat=subplot(1, 2, 1); 
imagesc(corrcoef(table2array(trainf(:, varidx)))), colormap(corrcmap'), axis equal tight, colorbar
% testmat=subplot(1, 2, 2); imagesc(corrcoef(table2array(test(:, varidx)))), colormap(corrcmap'), axis equal tight, colorbar
set(gca, 'clim', [-1 1], 'ytick', 1:length(varidx), 'yticklabel', lbls, 'xtick', 1:length(varidx), 'xticklabel', lbls, 'xticklabelrotation', 90, 'xaxislocation', 'top');
% set(testmat,  'clim', [-1 1], 'ytick', 1:length(varidx), 'yticklabel', lbls, 'xtick', 1:length(varidx), 'xticklabel', lbls, 'xticklabelrotation', 90, 'xaxislocation', 'top');

%% train madrs12
clc, tic, disp('Working...')
varidx=[4 9:23 25 26 28:30];
% varidx=[ 9:23 34 35];
for idx=1:10
    mdl=stepwiselm(trainf, 'ResponseVar', 'madrs12', 'PredictorVars', varidx, 'Criterion', 'aic', 'NSteps', idx, 'Lower', 'constant', 'Upper', 'interactions', 'Verbose', 0);
    % collect info on trained model
    mdlh(idx).varsin=mdl.PredictorNames; 
    mdlh(idx).weights=mdl.Coefficients; 
    mdlh(idx).performance=[mdl.Rsquared.Adjusted mdl.RMSE];
    
    % test performance on ket dataset
    mdlh(idx).test=[test.madrs12 feval(mdl, test(:, varidx))];
    r=corrcoef(test.madrs12, feval(mdl, test(:, varidx)));
    mdlh(idx).testp=[r(2)^2 sqrt(mean((test.madrs12-feval(mdl, test(:, varidx))).^2))];
end
toc

a=reshape([mdlh.performance], 2, [])';
b=reshape([mdlh.testp], 2, [])';
figure, plot(a(:, 2), a(:, 1), 'bo', b(:, 2), b(:, 1), 'ro');
legend({'trained', 'tested'});
set(gca, 'ylim', [0 1]);

%% plot correlation matrices madrs12
varidx=[6 4 9:23 25 26 28:30]; % select predictors
lbls=trainf.Properties.VariableNames(varidx); % get predictor names

lim=1; ramp=100;
corrcmap=[zeros(1, lim) linspace(0, 1, ramp) ones(1, 256-lim-ramp);...
    zeros(1, lim) linspace(0, 1, ramp) ones(1, 256-2*lim-2*ramp) linspace(1, 0, ramp) zeros(1, lim);...
    ones(1, 256-lim-ramp) linspace(1, 0, ramp) zeros(1, lim)];

figure, 
imagesc(corrcoef(table2array(trainf(:, varidx)))), colormap(gca, corrcmap'), axis equal tight
set(gca, 'clim', [-1 1], 'ytick', 1:length(varidx), 'yticklabel', lbls, 'xtick', 1:length(varidx), 'xticklabel', lbls, 'xticklabelrotation', 90, 'xaxislocation', 'top');


%% plot matrix of predictors 
c=[linspace(0, 1, 128) ones(1, 128); linspace(0, 1, 128) linspace(1, 0, 128); ones(1, 128) linspace(1, 0, 128)]; % make colormap blue-gray-red
[~, si]=sort(coefs1.RMSE, 'descend'); si=si([6:18 19 1:5]); % sort models by performance, then make look like Fig 2C
% varidx=[8 14 12 10 15 18 11  16 17 9]; % specify variables to show
varidx=[8 14 12 10 11 15 18 16 17]; % specify variables to show
data=table2array(coefs1(si, varidx)); % select data
subplot(1, 2, 1); imagesc(data), colormap(c'), colorbar, axis equal tight; % display heatmap
% make beautiful
set(gca, 'clim', [-2 2], 'xtick', 1:length(varidx), 'xticklabels', coefs1.Properties.VariableNames(varidx), 'XTickLabelRotation', 90, 'ytick', []); 

% add heatmap for intercept
subplot(1, 2, 2); imagesc(coefs1.intercept(si)), colormap(c'), colorbar, axis equal tight;
set(gca, 'clim', [-150 150], 'ytick', [], 'xtick', []);

%% try differences due to bin size in IS and IV
tic
rs=[1; unique(cumprod(perms(factor(1440)),2))]; % bin sizes for resampling
op=struct; % prepare output structure

for idx=1:length(cbt)
    sp=find(cbt(idx).tdata(:, 1)==0, 1, 'first');
    d=cbt(idx).tdata(sp:end, 2);
    pp=idvar(d, 1440, rs);
    op(idx).id=cbt(idx).ID; 
    op(idx).iv=pp.iv(:, 2)'; 
    op(idx).is=pp.is(:, 2)'; 
    clear pp sp d 
end
toc 

%% assess performance in dummy models 
a=[29;30;37;24;32;28;27;24;25;28;35;25]; % observed scores
b=repmat(mean(a), length(a), 1); % constant average
tic
n=[b min(a)+randi(14, [12, 1000000])-1]; % random numbers to match the range in observed scores

r=[sqrt(mean((a-n).^2)); mean(a-n); std(a-n)]; %[RMSE; mean error; sd error]
p=[sum(abs(a-n)<1); sum(abs(a-n)<2); sum(abs(a-n)<3); sum(abs(a-n)<4); sum(abs(a-n)<5)]/12; % [%<1, %<2, ...%<5]
toc

%% train madrs1 on all, test on test/train separately
clc, tic, disp('Working...')

varidx=[ 4 7 9 11:21];
for idx=1:12
    mdl=stepwiselm(trainf, 'ResponseVar', 'madrs1', 'PredictorVars', varidx, 'Criterion', 'aic', 'NSteps', idx, 'Lower', 'constant', 'Upper', 'interactions', 'Verbose', 0);
    % collect info on trained model
    mdlh(idx).varsin=mdl.PredictorNames; 
    mdlh(idx).weights=mdl.Coefficients; 
    mdlh(idx).performance=[mdl.Rsquared.Adjusted mdl.RMSE];
    r0=corrcoef(table2array(train(:, mdl.PredictorNames)));
    mdlh(idx).vif=diag(inv(r0)); % VIF in training dataset
    
    % test performance on ket dataset
    mdlh(idx).test=[trainf.madrs1 feval(mdl, trainf(:, varidx)) ];
    r=corrcoef(trainf.madrs1, feval(mdl, trainf(:, varidx)));
    mdlh(idx).testp=[r(2)^2 sqrt(mean((trainf.madrs1-feval(mdl, trainf(:, varidx))).^2))];
end
toc

a=reshape([mdlh.performance], 2, [])';
b=reshape([mdlh.testp], 2, [])';
figure, plot(a(:, 2), a(:, 1), 'bo'); %, b(:, 2), b(:, 1), 'ro');
% legend({'trained', 'tested'});
set(gca, 'ylim', [0 1]); xlim([0 8])
ylabel 'adjRsquare'; xlabel 'RMSE'

%% train madrs1 on depresjon madrs1>19
clc, tic, disp('Working...')

varidx=[ 4 7:21];
for idx=1:7
    mdl=stepwiselm(t, 'ResponseVar', 'madrs1', 'PredictorVars', varidx, 'Criterion', 'aic', 'NSteps', idx, 'Lower', 'constant', 'Upper', 'interactions', 'Verbose', 0);
    % collect info on trained model
    mdlh(idx).varsin=mdl.PredictorNames; 
    mdlh(idx).weights=mdl.Coefficients; 
    mdlh(idx).performance=[mdl.Rsquared.Adjusted mdl.RMSE];
    r0=corrcoef(table2array(t(:, mdl.PredictorNames)));
    mdlh(idx).vif=diag(inv(r0)); % VIF in training dataset
    mdlh(idx).train=[t.madrs1 mdl.Fitted];
    % test performance on cbt dataset
    mdlh(idx).test=[train.madrs1 feval(mdl, train(:, varidx)) ];
    r=corrcoef(train.madrs1, feval(mdl, train(:, varidx)));
    mdlh(idx).testp=[r(2)^2 sqrt(mean((train.madrs1-feval(mdl, train(:, varidx))).^2))];
end
toc

a=reshape([mdlh.performance], 2, [])';
b=reshape([mdlh.testp], 2, [])';
figure, plot(a(:, 2), a(:, 1), 'bo', b(:, 2), b(:, 1), 'ro');
legend({'trained', 'tested'});
set(gca, 'ylim', [0 1]); xlim([0 8])
ylabel 'adjRsquare'; xlabel 'RMSE'

%% brute force approach: train and test multiple regression model with restricted number of predictors
clc, tic, disp('Working...')
mdlh=struct('varsin', [], 'weights', [], 'train', [], 'performance', [], 'vif', [], 'test', [], 'testp', [], 'tvif', [], 'precision', []); % initialize output structure

varidx=[4 7:21]; % set of vars to consider for modelling; includes age and circadian period

for nidx=1:6 % limit to 6 predictors
varsel=nchoosek(1:length(varidx), nidx); % all possible combination of n vars
prg=waitbar(0, 'Working...');

for idx=1:length(varsel)
    waitbar(idx/length(varsel), prg);
    mdl=fitlm(train, 'ResponseVar', 'madrs1', 'PredictorVars', varidx(varsel(idx, :)));
    
    % collect info on trained model
    mdlh(end+1).varsin=mdl.PredictorNames; 
    mdlh(end).weights=mdl.Coefficients;
    mdlh(end).train=[train.madrs1 mdl.Fitted]; 
    mdlh(end).performance=[mdl.Rsquared.Adjusted mdl.RMSE];
    r0=corrcoef(table2array(train(:, mdl.PredictorNames)));
    mdlh(end).vif=diag(inv(r0)); % VIF in training dataset
    
    % test performance on test dataset
    [ypred, yci]=predict(mdl, test(:, varidx));
    mdlh(end).test=[test.madrs1 ypred yci];
    r=corrcoef(test.madrs1, ypred);
    r0=corrcoef(table2array(test(:, mdl.PredictorNames)));
    mdlh(end).testp=[r(2)^2 sqrt(mean((test.madrs1-ypred).^2))];
    mdlh(end).tvif=diag(inv(r0)); % VIF in test dataset
    mdlh(end).precision=[sum(abs(diff(mdlh(end).test(:, 1:2), 1, 2))<1); sum(abs(diff(mdlh(end).test(:, 1:2), 1, 2))<2);...
        sum(abs(diff(mdlh(end).test(:, 1:2), 1, 2))<3); sum(abs(diff(mdlh(end).test(:, 1:2), 1, 2))<4);...
        sum(abs(diff(mdlh(end).test(:, 1:2), 1, 2))<5)]/length(test.madrs1);
    mdlh(end).auc=sum(mdlh(end).precision);
    mdlh(end).slope=polyfit(test.madrs1, ypred, 1);
    
end % train next model
delete(prg); % kill progress bar
end % next 
toc

% crop out first line
mdlh=mdlh(2:end);
for idx=1:length(mdlh)
    mvif(idx)=max(mdlh(idx).vif); % get max VIF from each model
end
a=reshape([mdlh.performance], 2, [])'; % [Rsq RMSE] train
b=reshape([mdlh.testp], 2, [])'; % [Rsq RMSE] test

% identify models with acceptable accuracy RMSE<3 and precision
% R-square>0.5 and VIF<5 in train dataset; and RMSE<3 and R-square>.576 in test
% dataset

sidx=all([mvif(:)<5 a(:, 2)<=3 a(:, 1)>=.5 b(:, 2)<=3 b(:, 1)>=.576^2], 2); % find the models
[fv, ~, ia]=unique([cat(1, mdlh(sidx).varsin)]); % get variables involved

% put in table and sort by occurrence
t=table(fv, accumarray(ia, 1), 'variablenames', {'predictor', 'occurrence'});
t=sortrows(t, 'occurrence', 'descend'); 
t.occurrence=t.occurrence/sum(sidx); % show as percentage occurrence across models

% clean up 
clear fv ia idx nidx mdl prg r r0 varidx varsel yci ypred 
