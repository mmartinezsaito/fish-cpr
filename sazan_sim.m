function nll_a = sazan_sim(T, s, modtype, params)

% params: nummodels x numparams

datadir = '/home/mario/Yandex.Disk/CPR/logs/sim';

% utils
boltzfun = @(b, x, i) exp(-x(:, i) .* b) ./ sum(exp(-x .* b), 2);

% initialize variables
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
[nr nc] = size(T);

% subject-wise parameters
choice_t0 = hist(T{any(T.sid == s, 2), 'demanded'}, 3) / sum(T{any(T.sid == s, 2), 'demanded'});

nll = cell(ng, 1); % nll at each game
nll_a = 0; % nll summed over all trials
no = 0; % number of observations
R = cell(ng, nt); % generic reinforcement variable
E = R; % error values 
ME = R; % mean error values 


for i = 1:ng
  % initialize 
  fish = fish_t0; % amount of fish
  Q = 4 * choice_t0 .* ones(nm, 3); % Q value 
  ll_i = 0;  % log-likelihood
  
  % loop over all trials in the game
  for t = 1:nt  
    
    %demand = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'demanded'}; 
    %payoff = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'reaped'}; 
    shr1 = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'shrink1'}; 
    shr2 = T{any(T.sid == s, 2) & T.game == i & T.trial == t, 'shrink2'}; 
    fish = T{any(T.sid==s,2)&T.game==i&T.trial==t-1,'stock'};
    
    % calculate probabilities of choices and make a choice
    Pch = boltzfun(beta, -Q, net_sizes); 
    demand = randsample(3, 1, true, Pch);
    left_stock = fish - (shr1+shr2); 
    if     left_stock >= demand, payoff = demand;
    elseif left_stock >= 0,      payoff = left_stock;
    else,                        payoff = 0;      
    end
            
    % update the log-likelihood    
    ll_i = ll_i + log(Pch(:, demand)); 
    
    % RPE algorithms
    switch modtype
    case 'social'        
      comparison = net_sizes - mean([shr1 shr2]);
      R{i,t} = theta * payoff + (1-theta) .* comparison;
    case {'sustainable' 'nonsocial'}
      if size(params, 2) == 4 
        susmod2 = 1; kappa = params(:, 4); % ADD A NEW PARAMETER SO THAT sust CAn BE POSITIVE
      else
        susmod2 = 0 ; kappa = 0;
      end
      sust = kappa - abs(sustainable_popul - sum([shr1 shr2]) - net_sizes); 
      R{i,t} = theta * payoff + (1-theta) .* sust;        
    case 'modfreeRL'
      R{i,t} = repmat(payoff, nm, 1);         
    case 'inequity_aversion' 
      dia = params(:, 4); aia = params(:, 5);
      comparison = net_sizes - dia * sum(max([shr1 shr2] - payoff, 0)) ...
                             - aia * sum(max(payoff - [shr1 shr2], 0)); 
      R{i,t} = theta * payoff + (1-theta) .* comparison;        
    end
      
    rpe = R{i,t} - Q; 
    rpe(setdiff(net_sizes, demand)) = 0; 
      
    if wrtvars
      T{any(T.sid == s, 2) & T.game == i & T.trial == t, ['rpe_sim' modtype(1:6)]} = rpe;
      T{any(T.sid == s, 2) & T.game == i & T.trial == t, ['val_sim' modtype(1:6)]} = Q;
      T{any(T.sid == s, 2) & T.game == i & T.trial == t, ['pch_sim' modtype(1:6)]} = Pch;
    end
      
    Q = Q + alpha .* rpe;    
               
  end % j
  
  nll{i} = -ll_i;
  nll_a = -ll_i + nll_a; 
end % i


switch modtype
case 'inequity_aversion',   fprintf(1, '%u %f %f %f %f %f %f\n', s, nll_a, params)
case 'modfreeRL',           fprintf(1, '%u %f %f %f\n', s, nll_a, params(1:2))
case 'nonsocial' & susmod2, fprintf(1, '%u %f %f %f %f %f\n', s, nll_a, params)
otherwise,                  fprintf(1, '%u %f %f %f %f\n', s, nll_a, params)
end



