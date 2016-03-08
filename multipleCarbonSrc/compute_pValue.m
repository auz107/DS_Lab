clear all

% Experimental observations (1 = A dominates,  2 = B dominates, 3 = A and B coexist)
exp = [1,3,3,3,3,3,3,3,3,2,3,2,3,3,2];

% A matrix of 10000 random permutations of exp
rand_mat = zeros(10000,15);
for i = 1:10000
    rand_mat(i,:) = exp(randperm(15));
end

rand_mat = randi([1,3],[10000,15]);

% Count the number of matches
for i = 1:10000
    match_nums(i) = length(find(rand_mat(i,:) == exp));
end

% Plot the histogram
%hist(match_nums)

% Mean and std 
mean_val = mean(match_nums)
std_val = std(match_nums)

% z-score
[Z,mu,sigma] = zscore(match_nums);
mean(normcdf(-abs(Z),mu,sigma))

% Normal cummulative distirbution 
%cdf = normcdf(match_nums,mean_val,std_val);  

% Run test
[h,p] = runstest(match_nums,7)  
