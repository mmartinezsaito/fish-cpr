function P = sazan_optim(T, ot)

soc_sids = [2:10 20 21 26 27 30 31 34 35 37 38 44 45 47 49 51];
nsoc_sids = setdiff(2:51, soc_sids);

bic = @(n, k, nll) k * log(n) + 2 * nll;

ss = unique(T.sid)';
ns = length(ss);
xin = [0.5 1 0.5];

modtypes = {'social' 'nonsocial' 'modfreeRL' 'inequity_aversion'};

fmcot = {'interior-point', 'sqp', 'active-set', 'interior-point', 'interior-point'}; % 'trust-region-reflective'
gs = GlobalSearch;
ms = MultiStart;
ds_options = optimoptions('patternsearch');
ps_options = optimoptions('particleswarm');
sa_options = saoptimset('simulannealbnd');
ga_options = gaoptimset('MutationFcn', @mutationadaptfeasible);

p = {}; P = [];
for m = 0:4  %<--- RESTORED
  switch m
  case 4,    x0 = [xin 0 0]; lb = [0 0 0 0 0];  ub = [1 inf 1 1 1];
  %case 2,    x0 = [xin 0];   lb = [0 0 0 -inf]; ub = [1 inf 1 inf];
  otherwise, x0 = xin;       lb = [0 0 0];      ub = [1 inf 1];
  end  
  
  for s = ss      
    no = length(T{T.sid == s, 'demanded'}); 
    if m ~= 0  
      fprintf(1, 'subject\tnull\tm_nll\talpha\tbeta\t(theta)\t(dia)\t(aia)\n')  
      
      switch ot
      case {1 2 3}
        fmc_options = optimoptions('fmincon', 'Algorithm', fmcot{ot});
        [p, nll, exflag(s, m), output(s, m), lambda_p, grad_p, hess_p] = ...
          fmincon(@(x) sazan_nll(T, s, modtypes{m}, x, 0), ...
          x0, [], [], [], [], lb, ub, [], fmc_options)
      case 4
        ops = createOptimProblem('fmincon', 'objective', ...
          @(x) sazan_nll(T, s, modtypes{m}, x, 0), 'x0', x0, 'lb', lb, 'ub', ub, ...
          'options', optimoptions('fmincon', 'Algorithm', fmcot{ot}));
        [p, nll, exflag(s, m), output(s, m), sol] = run(gs, ops);
      case 5
        ops = createOptimProblem('fmincon', 'objective', ...
          @(x) sazan_nll(T, s, modtypes{m}, x, 0), 'x0', x0, 'lb', lb, 'ub', ub, ...
          'options', optimoptions('fmincon', 'Algorithm', fmcot{ot}));
        [p, nll, exflag(s, m), output(s, m), sol] = run(ms, ops, 10);
      case 6
        [p, nll, exflag(s, m), output(s, m)] = patternsearch( ...
          @(x) sazan_nll(T, s, modtypes{m}, x, 0), ...
          x0, [], [], [], [], lb, ub, [], ds_options);
      case 7
        [p, nll, exflag(s, m), output(s, m)] = particleswarm( ...
          @(x) sazan_nll(T, s, modtypes{m}, x, 0), length(x0), lb, ub, ps_options);
      case 8
        [p, nll, exflag(s, m), output(s, m)] = simulannealbnd( ...
          @(x) sazan_nll(T, s, modtypes{m}, x, 0), x0, lb, ub, sa_options);          
      case 9
        [p, nll, exflag(s, m), output(s, m), popl, scores] = ga( ...
            @(x) sazan_nll(T, s, modtypes{m}, x, 0), length(x0), ...
            [], [], [], [], lb, ub, [], ga_options);          
      end
    
    else
      nll =  no * log(3);
    end
      
    switch m
    case 4,     params = p;                    k = 5;
    case 3,     params = [p(1:2) NaN NaN NaN]; k = 2; 
    case 2,     params = [p NaN NaN];          k = 3;  % [p NaN]; k = 4; 
    case 1,     params = [p NaN NaN];          k = 3;  
    case 0,     params = [1 NaN NaN NaN NaN];  k = 1;
    end
    
    P = [P; table(m, s, no, nll, bic(no, k, nll), params, ...
        'VariableNames', {'model' 'sid' 'numobs' 'nll' 'bic' 'params'})];
  end
end
P.issoc = ismember(P.sid, soc_sids);
P.model = categorical(P.model, 'Ordinal', 0);

