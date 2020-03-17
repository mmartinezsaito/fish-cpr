if isunix,   bdir = '/home/mario/Data/Yandex.Disk/CPR/'; 
elseif ispc, bdir = 'C:\Users\Administrator\Yandex.Disk\CPR\'; 
end
datadir = [bdir 'behdat'];

start_time = '1stpulse'; % '1stpulse' '1stsession' 
TR = 2.28; %s

soc_sids = [2:10 20 21 26 27 30 31 34 35 37 38 44 45 47 49 51];
nsoc_sids = [11:19 22:25 28 29 32 33 36 39 40:43 46 48 50];

cd(datadir)
ddir = dir(datadir);
k = 0; fn = {};
for i = 1:length(ddir)
  if strfind(ddir(i).name, 'testing')
    k = k + 1;  
    fn{k} = ddir(i).name; %fn{k}.id = fopen([datadir filesep fn{k}.name], 'r');        
  end
end

formatSpec = '%*s %u %s %s %f %*u %*u %*u %*u %*u %*s';

T0 = {}; T1 = {}; prev_tmark = 0;
for j = 1:k    
  T0{j} = readtable(fn{j}, 'FileType', 'text', 'HeaderLines', 5, ...
    'Format', formatSpec, 'Delimiter', '\t', 'ReadVariableNames', 0);
  T0{j}.Properties.VariableNames = {'trial', 'event', 'code', 'time'};
  firstpulsetime = T0{j}.time(1) / 10^4;
  T0{j} = T0{j}(~strcmp(T0{j}.event, 'Pulse'), :);
  
  nt = size(T0{j},1);
  sid = regexp(fn{j}, '^s(.*?)_', 'tokens'); sid = str2num(sid{1}{1});
  T0{j}.sid = sid * ones(nt, 1);
  
  switch start_time
  case '1stsession' 
    for i = 1:nt        
      if strcmp(T0{j}.code(i), 'New Session')  
        synct(j) = T0{1}.time(i) / 10^4; break
      end    
    end
  case '1stpulse'
    synct(j) = firstpulsetime + 4 * TR;  % +4TR because vasily did so
  end  
  T0{j}.time = T0{j}.time / 10^4 - synct(j);
     
  T0{j}.game = NaN(nt, 1);     T0{j}.trial = NaN(nt, 1); 
  T0{j}.demanded = NaN(nt, 1); T0{j}.reaped = NaN(nt, 1); 
  T0{j}.shrink1 = NaN(nt, 1);  T0{j}.shrink2 = NaN(nt, 1); T0{j}.stock = NaN(nt, 1); 
  T0{j}.shr_time = NaN(nt, 1); T0{j}.dem_time = NaN(nt, 1); T0{j}.end_time = NaN(nt, 1);
  
  pmark = 1; 
  for i = 1:nt
    S = regexp(T0{j}.code{i}, '^([0-9]){1,2}\(([0-9])\)$', 'tokens');
    if ~isempty(S)
      if pmark ~= 1
        T0{j}.trial(pmark:i-1) = tmark*ones(i-pmark,1);
        T0{j}.game(pmark:i-1) = smark*ones(i-pmark,1);
      end
      pmark = i+1; tmark = str2double(S{1}{2}); smark = str2double(S{1}{1});   
    elseif ~isempty(regexp(T0{j}.code{i}, '^End', 'once')) && smark == 16
      T0{j}.trial(pmark:i-1) = tmark*ones(i-pmark,1);
      T0{j}.game(pmark:i-1) = smark*ones(i-pmark,1);
    end
    
    Sc = regexp(T0{j}.code{i}, ...
        '^(\w+): ([0-9])$|^(\d) (\d)$|^([0-9])$|(F)(\d{1,2})|^(End), .*', 'tokens');
    if ~isempty(Sc)
      if prev_tmark ~= tmark || fi > size(T0{j},1), fi = i; prev_tmark = tmark; end
    
      if strcmp(T0{j}.event{i}, 'Response') && i > 2
        T0{j}.dem_time(fi) = T0{j}.time(fi-2); continue       
      end
      
      switch Sc{1}{1}
      case 'Own'
        T0{j}.reaped(fi) = str2double(Sc{1}{2});
        T0{j}.shr_time(fi) = T0{j}.time(fi+2); % IMPORTANT: shr_time is feedback about shrink event 
      case 'Demanded'
        T0{j}.demanded(fi) = str2double(Sc{1}{2}); 
      case {'0' '1' '2' '3'}
        if length(Sc{1})~=2, continue; end
        T0{j}.shrink1(fi) = str2double(Sc{1}{1}); 
        T0{j}.shrink2(fi) = str2double(Sc{1}{2});   
      case {'F' 'End'}
        T0{j}.end_time(fi) = T0{j}.time(fi+5);              
        if length(Sc{1}) == 2
          T0{j}.stock(fi) = str2double(Sc{1}{2});              
        else
          T0{j}.stock(fi) = 0;              
        end
      end
      
    end
  end  
  
  T1 = [T1; T0{j}(~isnan(T0{j}.reaped), [1 4:end])]; % time is for reaped fish feedback
end
T1.issoc = ismember(T1.sid, soc_sids);
T1.Properties.VariableNames{2} = 'rea_time';
T1 = T1(:, [3 4 1 5:9 11 2 10 12]);

T = T1;

cd(bdir)

