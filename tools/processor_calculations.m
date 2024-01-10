%% Script to decide number of nodes and or domain size

X = 1152;
Y = 768;
Z = 384;
ppc_lim = 360000;
nodes = 8;
cpernode = 128;
tot_points = X*Y*Z;
proxs = alldivisors(X);
proys = alldivisors(Y);
tot_pros = proxs.*proys';
pairs = [];
min_nodes = [];
for i = 1:length(proxs)
    for j = 1:length(proys)
        prox = proxs(i);
        proy = proys(j);
        ppc = tot_points/(prox*proy);
        if ppc < ppc_lim
            %if ((X/prox)*(Y/proy)) > 54
                minnode = (prox*proy)/cpernode;
                min_nodes = [min_nodes, minnode];
                if prox*proy <= nodes*cpernode
                     pairs = [pairs; prox,proy];
                end 
            %end 
        end 
    end 
end
function divs = alldivisors(N)
  % compute the set of all integer divisors of the positive integer N
  
  % first, get the list of prime factors of N. 
  facs = factor(N);
  
  divs = [1,facs(1)];
  for fi = facs(2:end)
    % if N is prime, then facs had only one element,
    % and this loop will not execute at all, In that case
    % The set of all divisors is simply 1 and N.
    
    % this outer product will generate all combinations of
    % the divisors found so far, combined with the current
    % divisor fi.
    divs = [1;fi]*divs;
    
    % unique eliminates the replicate divisors, making
    % this an efficient code.
    divs = unique(divs(:)');
  end
  
end