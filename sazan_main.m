%% read logs
readlog_sazan


%% optimize using 9 different solvers 
otn = {'fmc_ip', 'fmc_sqp', 'fmc_as', 'gs_fmc_ip', 'ms_fmc_ip', ...
      'patsrch', 'pswarm', 'sabnd', 'ga'}; 
for ot = [1:4 6:9]  % skip 5
  Pt = sazan_optim(T, ot);
  if ot==1
    P = Pt; 
    P.Properties.VariableNames{'nll'} = ['nll_' otn{ot}]; 
    P.Properties.VariableNames{'bic'} = ['bic_' otn{ot}]; 
    P.Properties.VariableNames{'params'} = ['pars_' otn{ot}]; 
  else
    P{:, ['nll_' otn{ot}]} = Pt.nll;
    P{:, ['bic_' otn{ot}]} = Pt.bic;
    P{:, ['pars_' otn{ot}]} = Pt.params;
  end
end
%writetable(P,'P.xlsx', 'FileType', 'spreadsheet', 'WriteVariableNames',true)

%% find best parameters
bic_cols = [5 9 12 15 18 21 24 27]%30];
P0 = P(P.model~=categorical(0), :);
nLmat = P{:, bic_cols};
[~, I] = sort(nLmat, 2);
tabulate(I(:,1))

sum(nLmat == min(nLmat, [], 2))
[idx jdx] = find(I' == 1);
P.bic_idx = idx;
P.bic_min = zeros(250,1); P.pars_best = zeros(250,5);
for i = 1:250
  P.bic_min(i) = P{i, bic_cols(idx(i))};
  P.pars_best(i,:) = P{i, bic_cols(idx(i))+1};
end
P.bic_min_pern = P.bic_min ./ P.numobs;

%P(P.model~=categorical(0),[1:3 end-2:end])

%% write dm models proxy vars using fitted parameters
modtypes = {'social' 'nonsocial' 'modfreeRL' 'inequity_aversion'};

for m = categorical(1:4)
  for s = unique(T.sid)'
    T = sazan_nll(T, s, modtypes{m}, P{P.model == m & P.sid == s, 'pars_best'}, 1); 
  end
end
%writetable(T,'T.xlsx', 'FileType', 'spreadsheet', 'WriteVariableNames',true)

%% plot

subplot(131), boxplot(P.nll_fmc_ip, P.model, 'jitter', 0.1, 'notch', 'on', 'datalim', [-inf inf], 'extrememode', 'compress', 'whisker', 1.5)
subplot(132), boxplot(P.nll_fmc_as, P.model, 'jitter', 0.1, 'notch', 'on', 'datalim', [-inf inf], 'extrememode', 'compress', 'whisker', 1.5)
subplot(133), boxplot(P.nll_fmc_sqp, P.model, 'jitter', 0.1, 'notch', 'on', 'datalim', [-inf inf], 'extrememode', 'compress', 'whisker', 1.5)

% BIC by sid
boxplot(P.bic_min_pern, P.sid - 1, 'colorgroup', P.issoc, 'plotstyle', 'traditional')
xlabel('Subject ID'), ylabel('BIC score per observation'), ylim([.6 2.6]) %title('BICs from the four models') 

% fig03: BIC by model and issoc
boxplot(P.bic_min_pern, {P.issoc P.model}, 'ColorGroup', P.issoc, ...
    'Labels', {num2cell('nnnnnsssss'); num2cell('0snri0snri')}, ...  %{mat2cell(['nnnnnsssss'; '0snri0snri'],[1 1],ones(1,10))}
    'jitter', 0.1) %, 'notch', 'on')
xlabel('sociality treatment and model type'), ylabel('BIC score per observation'), %title('Goodness of fit by model type and treatment') 

P2 = P(P.model==categorical(1) | P.model==categorical(2), :);
boxplot(P2.bic_min_pern, {P2.issoc P2.model})
xlabel('treatment and model type'), ylabel('BIC score per observation'), title('Goodness of fit by model type and treatment') 

% figS2: BIC difference between nsoc and soc models by issoc
boxplot(P{P.model==categorical(2), 'bic_min_pern'} - P{P.model==categorical(1), 'bic_min_pern'}, P{P.model==categorical(1),'issoc'}, ...
    'Labels', {num2cell('ns')}) %, 'notch', 'on')
xlabel('condition'), ylabel('BIC score difference (n-s) per observation'), %title('BIC scores differences between nonsocial and social treatments') 

% figS3: parameters
P3 = P(P.model~=categorical(0),:);
pcell={'alpha','beta','theta','dia','aia'}; figure
for i = 1:5
  subplot(1,5,i) %subplot(3,5,i) 
  boxplot(P3.pars_best(:,i), {P3.issoc P3.model}, 'jitter', 0.1, 'datalim', [0 20], ...
      'extrememode', 'compress', 'whisker', 1.5, 'labels', {num2cell('nnnnssss'); num2cell('snrisnri')}) %,'notch', 'on') 
  title(pcell{i})
  %subplot(3,5,i+5) 
  %histogram2(int8(P(P.issoc,:).model), P(P.issoc,:).pars_best(:,i), [10 25], 'FaceAlpha', 0.3)  
  %subplot(3,5,i+10) 
  %histogram2(int8(P(~P.issoc,:).model), P(~P.issoc,:).pars_best(:,i), [10 25], 'FaceAlpha', 0.3)  
end


%% simple stats

% average BIC per observation
[m, s]=grpstats(P.bic_min_pern, P.model, {'mean', 'std'})


% between subjects num of trials per condition
soc_sids = [2:10 20 21 26 27 30 31 34 35 37 38 44 45 47 49 51];
nsoc_sids = [11:19 22:25 28 29 32 33 36 39 40:43 46 48 50];

issoc = ismember(T.sid, soc_sids);
sT=T(issoc,:);
nT=T(~issoc,:);
nts = grpstats(sT.trial, sT.sid, {'numel'}) / 16;
ntn = grpstats(nT.trial, nT.sid, {'numel'}) / 16;
mean(nts), std(nts)
mean(ntn), std(ntn)
[H,P,CI,STATS] = ttest2(nts,ntn)



%% fit linear models

% linear mixed effects model
lme = fitlme(P, 'bic_min_pern ~ model*issoc + (1|sid)', ...
      'FitMethod', 'REML', 'Optimizer', 'quasinewton', ...
      'CheckHessian', true, 'CovariancePattern', 'FullCholesky', ...
      'DummyVarCoding', 'reference');
lme                       % Table 2
plotResiduals(lme)
%coefCI(lme)
%randomEffects(lme)'
anova(lme, 'DFMethod', 'satterthwaite') %residual 
%lme.Rsquared
lme.CoefficientCovariance


% fixed effects anova: not significant
[pv, t, stats, terms] = anovan(P2.bic_min_pern, [int8(P2.model) P2.issoc], 'model', 2, 'varnames', {'mty', 'issoc'});
[comp, mns, h, gnames] = multcompare(stats, 'alpha', .05, 'ctype', 'tukey-kramer')


% mixed anova of soc and nsoc respect to null 
Ps0 = P(ismember(uint8(P.model), [1 2]), :);
Pn0 = P(ismember(uint8(P.model), [1 3]), :);
[P,T,STATS,TERMS] = anovan(Ps0.bic_min_pern, {Ps0.model, Ps0.sid}, 'random', 2)
[P,T,STATS,TERMS] = anovan(Pn0.bic_min_pern, {Pn0.model, Pn0.sid}, 'random', 2)


% Mixed anova
brmP = table(P{1:50,32}, P{51:100,32}, P{101:150,32}, P{151:200,32}, P{201:250,32}, P{1:50,7}, 'VariableNames', {'z','s','n','r','i', 'issoc'});
wrmP = table([0:4]', 'VariableNames', {'model'});
rm = fitrm(brmP, 'z,s,n,r,i ~ issoc', 'WithinDesign', wrmP, 'WithinModel', 'separatemeans');
brmP2 = table(P{51:100,32}, P{101:150,32}, P{1:50,7},'VariableNames', {'s','n', 'issoc'});
wrmP2 = table([1 2]', 'VariableNames', {'model'});
rm2 = fitrm(brmP2, 's,n ~ issoc', 'WithinDesign', wrmP2, 'WithinModel', 'separatemeans');

anovatbl = anova(rm) % anova for between-subject effects
mauchly(rm) % sphericity test
epstbl = epsilon(rm)
ranovatbl = ranova(rm) % rmanova (for within-subject effects)
mctbl = multcompare(rm, 'model', 'alpha', .05, 'ComparisonType', 'tukey-kramer')

% to report: Table3
anovatbl = anova(rm2)
ranovatbl = ranova(rm2)