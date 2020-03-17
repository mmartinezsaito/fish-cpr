function nll_m = sazan_nll(T, s, modtype, params, wrtvars)

% params: nummodels x numparams

if isunix,   datadir = '/home/mario/Yandex.Disk/CPR/logs';
elseif ispc, datadir = 'C:\Users\Administrator\Yandex.Disk\CPR\logs';    
end

% utils
boltzfun = @(b, x, i) exp(-x(:, i) .* b) ./ sum(exp(-x .* b), 2);

% parameters
iscongruent = 1;

fish_t0 = 18;  % initial amount of fish in the lake
recov_rate = 1.5; % fish recovery rate
net_sizes = [1 2 3]; % possible choices 
sustainable_popul = 6;

alpha = params(:, 1); % learning        
beta = params(:, 2);  % inverse temperature
theta = params(:, 3); % fairness coefficient

np = size(params, 2);
nt = 8; % trials in one game
ng = 16; % games in one session
nm = size(params, 1); % number of models
[nr nc] = size(T);


% subject-wise parameters
%if s == 0, s = unique(T.sid)'; end

choice_t0 = hist(T{any(T.sid == s, 2), 'demanded'}, 3) / sum(T{any(T.sid == s, 2), 'demanded'});

nll = cell(ng, 1); % nll at each game
nll_m = zeros(nm, 1); % nll summed over all trials for each model
nll_null = 0; % baseline model
no = 0; % number of observations
R = cell(ng, nt); % generic reinforcement variable
E = R; % error values 
ME = R; % mean error values 


for i = 1:ng
  % initialize 
  fish = fish_t0; % amount of fish
  Q = 4 * choice_t0 .* ones(nm, 3); % Q value 
  ll_i = zeros(nm, 1);  % log-likelihood
    
  % loop over all trials in the game
  for t = 1:nt  
    demand = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'demanded'}; 
    payoff = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'reaped'}; 
    shr1 = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'shrink1'}; 
    shr2 = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'shrink2'}; 
      
    % calculate probabilities of choices
    Pch = boltzfun(beta, -Q, [1 2 3]); 
        
    % update the log-likelihood    
    if ismember(demand, 1:3) 
      ll_i = ll_i + log(Pch(:, demand)); 
      nll_null = nll_null - log(1/3);
      no = no + 1;
    end
       
    % update the amount of fish real in the lake
    fish_afterSeize = fish - demand - shr1 - shr2;
    if fish_afterSeize < 0, fish_afterSeize = 0; end
    fish = round(fish_afterSeize * recov_rate);
        
    if ismember(demand, 1:3) 
      switch modtype
      case 'social'
        comparison = net_sizes - mean([shr1 shr2]);
        R{i,t} = theta * payoff + (1-theta) .* comparison;
      case {'sustainable' 'nonsocial'}
        if size(params, 2) == 4 
            susmod2 = 1; kappa = params(:, 4); % ADD A NEW PARAMETER SO THAT sust CAN BE POSITIVE
        else
            susmod2 = 0 ; kappa = 0;
        end
        sust = kappa -abs(sustainable_popul - sum([shr1 shr2]) - net_sizes); 
        R{i,t} = theta * payoff + (1-theta) .* sust;        
      case 'modfreeRL'
        R{i,t} = repmat(payoff, nm, 1);         
      case 'inequity_aversion' 
        dia = params(:, 4); aia = params(:, 5);
        comparison = net_sizes - dia * sum(max([shr1 shr2] - payoff, 0)) ...
                               - aia * sum(max(payoff - [shr1 shr2], 0)); 
        R{i,t} = theta * payoff + (1-theta) .* comparison;        
      end
      
      % RPE computation
      rpe = R{i,t} - Q;                             % using counterfactual RPE
      %rpe = R{i,t}-Q; rpe(setdiff(net_sizes, demand)) = 0;         % local RPE 
      
      if wrtvars
        T{any(T.sid == s, 2) & T.game == i & T.trial == t, ['rpe_' modtype(1:6)]} = rpe;
        T{any(T.sid == s, 2) & T.game == i & T.trial == t, ['val_' modtype(1:6)]} = Q;
        T{any(T.sid == s, 2) & T.game == i & T.trial == t, ['pch_' modtype(1:6)]} = Pch;
      end
      
      Q = Q + alpha .* rpe;    %  Q(demand) = Q(demand) + alpha .* rpe(demand);    
    end   
        
  end % j
  
  nll{i} = -ll_i;
  nll_m = -ll_i + nll_m; 
end % i


%s
%no
%R, E, ME
%nll_null % log(3) * no
%nll_m
nll_best = min(nll_m);
params_best = params(nll_best == nll_m, :);
if length(params_best) ~= np, params_best = NaN(1, np); end

switch modtype
case 'inequity_aversion',     fprintf(1, '%u %.2f %f %f %f %f %f %f\n', s, nll_null, nll_m, params_best)
case 'modfreeRL',             fprintf(1, '%u %.2f %f %f %f\n', s, nll_null, nll_m, params_best(1:2))
case 'nonsocial', if susmod2, fprintf(1, '%u %.2f %f %f %f %f %f\n', s, nll_null, nll_m, params_best(1:4))
                  else,       fprintf(1, '%u %.2f %f %f %f %f %f\n', s, nll_null, nll_m, params_best(1:3)), end
otherwise,                    fprintf(1, '%u %.2f %f %f %f %f\n', s, nll_null, nll_m, params_best)
end

if wrtvars, nll_m = T; end





%{
old version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Stores only the final Log-Likelihood.                                  %%
%%Reads the desired participant data, predicts every next choice of the  %%
%%participant. Compares the predicted values with the real experimental  %%
%%data.                                                                  %%
%%Enables to modify only one free parameter of the model: theta           %%
%%Input Variables:                                                       %%
%%...@social:1=social condition, 0=nonsocial condition                   %%
%%...@sn: Participant number                                            %%
%%...@theta: freee parameter, fairness coefficient (range 0:1)            %%
%%...@alpha: freee parameter, learning rate (range 0:1)                  %%
%%...@beta: freee parameter, sensitivity parameter (range 1:100)        %%
%%Output Variable:                                                       %%
%%...@AllRoundsLL:  Log-likelihoods of the model summed over all models  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define which participant data to read
if iscongruent, datadir = [datadir '/processed_final'];
else,           datadir = [datadir '/processed_final_reverse'];
end

if issocial    
  filename = [datadir '/SOC/s' num2str(sn) '.txt'];
else
  filename = [datadir '/NSOC/s' num2str(sn) '.txt'];
end

%%Read the data from the desired datafile
data = load(filename);
data = sortrows(data, 1);

%% Initialise fixed variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of games, inferred separately for each participant
n = data(end, 1);
t = 8;

%% Initialise result variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Initialise result matrix
ModelLLs = cell(n,1); %Model log-likelihoods at each step
AllRoundsLL=zeros(size(params(:,1)));
Baseline = 0; %Baseline model
N=0; %Amount of observations
RMat=99.*ones(n,t); %Array with all R values per participant 
Error=RMat; %Array with error values per participant per trial
MeanError=RMat; %Array with mean error values per parthttps://www.youtube.com/watch?v=DDEQusMDhy8&index=766&list=PL01A9AA4DE2B6FEEAicipant per trial


%% Rearrange data from data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
competitorChoices = data(:,4:5); %this one has to be transposed
participantChoices = data(:,3);
participantReceived = data(:,2);
games=data(:,1);

%%Process the data file
%%Initialise new data matrices
pptChoices = zeros(n,t);
pptReceived = zeros(n,t);
compChoices1 = zeros(n,t);
compChoices2 = zeros(n,t);

for i=1:n
    index=(games==i);
    if (sum(index)>t)
        x=sum(index)-t;
        newIndex=find(index, x,'last');
        index(newIndex)=0;
    end
    pptChoices(i,:)=[participantChoices(index)', zeros(1,t-length(participantChoices(index)))];
    pptReceived(i,:)=[participantReceived(index)', zeros(1,t-length(participantReceived(index)))];
    compChoices1(i,:)=[competitorChoices(index,1)', zeros(1,t-length(competitorChoices(index,1)))];
    compChoices2(i,:)=[competitorChoices(index,2)', zeros(1,t-length(competitorChoices(index,2)))];
end


%% Find the most frequent starting values for the participant %%%%%%%%%%%%%
freq1=sum(pptChoices(:,1)==1)/n;
freq2=sum(pptChoices(:,1)==2)/n;
freq3=sum(pptChoices(:,1)==3)/n;

FREQS=[freq1, freq2, freq3];

%% Loop over all trials
%%Loop over all games
for i=1:n
    %%Initialise the amount of fish
    fishR = fish_t0; %Real fish
    
    %%Set the value before the first choice
    startQ=[4*freq1.*ones(size(params(:,1))), 4*freq2.*ones(size(params(:,1))), 4*freq3.*ones(size(params(:,1)))];
    
    %%Initialise the Model likelihood
    ModelLL = zeros(length(params),1);
    %%Loop over all trials in the game
    for j=1:t
        %%Initialise the data matrices
        prob = zeros(length(params),3); %probability of the choice
        Q = zeros(length(params),3); %Expected value of a choice
        
        %%Calcluate probabilities of choices
        prob(:,1)=1./(1+exp(beta.*(startQ(:,3)-startQ(:,1)))+exp(beta.*(startQ(:,2)-startQ(:,1))));
        prob(:,2)=1./(1+exp(beta.*(startQ(:,3)-startQ(:,2)))+exp(beta.*(startQ(:,1)-startQ(:,2))));
        prob(:,3) = 1 - sum(prob(:,1:2), 2);
        
        
        %%Update the log-likelihood, depending on the real data
        if pptChoices(i,j)==1
             ModelLL = ModelLL + log(prob(:,1));
             Baseline=Baseline + log(1/3);
             N=N+1;
        elseif pptChoices(i,j)==2
             ModelLL = ModelLL + log(prob(:,2));
             Baseline=Baseline + log(1/3);
             N=N+1;
        elseif pptChoices(i,j)==3
             ModelLL = ModelLL + log(prob(:,3));
             Baseline=Baseline + log(1/3);
             N=N+1;
        end
        
        
        %%Update the amount of fish real in the lake
        fishAfterR = fishR-pptChoices(i,j)-compChoices1(i,j)-compChoices2(i,j);
        if fishAfterR<=0
            fishAfterR = 0;
        end
        fishR = round(fishAfterR*recov_rate);
        if ismember(pptChoices(i,j), [1,2,3])
            Comparison=zeros(3,1);
            R=zeros(length(params),3);
            if issocial==1
                %Comparison=pptChoices(i,j)-mean([compChoices1(i,j), compChoices2(i,j)]);
                Comparison=(net_sizes-mean([compChoices1(i,j), compChoices2(i,j)]));
 
                R(:,1)=theta.*pptReceived(i,j)+(1-theta).*Comparison(1);
                R(:,2)=theta.*pptReceived(i,j)+(1-theta).*Comparison(2);
                R(:,3)=theta.*pptReceived(i,j)+(1-theta).*Comparison(3);
            elseif issocial==0
                %Comparison=-abs(6.*ones(1,3)-sum([compChoices1(i,j), compChoices2(i,j)])-choices); %LAST VALID SO FAR IN THE PAPER

                
                Comparison=-((6-sum([compChoices1(i,j), compChoices2(i,j)])-net_sizes).^3)+6; %Rsult 2 
                %Comparison=3.*(1./((6.0000001-sum([compChoices1(i,j), compChoices2(i,j)])-choices).^4))-1;
                %Comparison=(choices+sum([compChoices1(i,j), compChoices2(i,j)])-6);
                
                R(:,1)=theta.*pptReceived(i,j)+(1-theta).*Comparison(1);
                R(:,2)=theta.*pptReceived(i,j)+(1-theta).*Comparison(2);
                R(:,3)=theta.*pptReceived(i,j)+(1-theta).*Comparison(3);
            end            
            Q(:,1)=startQ(:,1)+alpha.*(R(:,1)-startQ(:,1));
            Q(:,2)=startQ(:,2)+alpha.*(R(:,2)-startQ(:,2));
            Q(:,3)=startQ(:,3)+alpha.*(R(:,3)-startQ(:,3));
        
            startQ = Q;
        end
    end
    ModelLLs{i}=-2.*ModelLL;
    AllRoundsLL=-2.*ModelLL+AllRoundsLL; %Sum of the LL
end

Baseline=-2*Baseline;
bestLL=min(AllRoundsLL);
indexBest=(bestLL==AllRoundsLL);
bestParams=(params(indexBest,:));
%}