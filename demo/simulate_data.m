function t = simulate_data()
% This function will simulate data that can be used to demo the SILA
% algorithm. A simple model that assumes slope = zero for non-accumulators
% and a fixed linear slope of 5 per year above the positivity threshold of
% 21.

%% setup the initial subjects
h=1;
while h~=0
    age_base_neg = 65+ 15*rand(130,1);
    age_base_pos = 65+ 15*rand(70,1);
    h = ttest2(age_base_neg,age_base_pos);
end

%% simluate data for negative cases
subid_neg = [1:130]';
mneg = 0; sdneg = 7;
vals_neg = normrnd(mneg,sdneg,numel(subid_neg)*3,1);

subid_neg_tall = nan(numel(subid_neg)*3,1);
subage_neg = nan(numel(subid_neg)*3,1);
subval_neg = nan(numel(subid_neg)*3,1);

for i = 1:numel(subid_neg)
    idx_low = 3*i-2;
    idx_hi = idx_low+2;
    subid_neg_tall(idx_low:idx_hi) = subid_neg(i);
    subage_neg(idx_low:idx_hi) = [age_base_neg(i);age_base_neg(i) + rand(2,1)-.5 + [2;4]];
    subval_neg(idx_low:idx_hi) = vals_neg(idx_low:idx_hi);
    
end

%% simulate data for positive cases
subid_pos = [131:200]';
val_base_pos = 100*rand(70,1)+5;
sd_pos = normrnd(0,sdneg,numel(subid_pos)*3,1);

subid_pos_tall = nan(numel(subid_pos)*3,1);
subage_pos = nan(numel(subid_pos)*3,1);
subval_pos = nan(numel(subid_pos)*3,1);

for i = 1:numel(subid_pos)
    idx_low = 3*i-2;
    idx_hi = idx_low+2;
    subid_pos_tall(idx_low:idx_hi) = subid_pos(i);
    subage_pos(idx_low:idx_hi) = [age_base_pos(i);age_base_pos(i) + rand(2,1)-.5 + [2;4]];
    % more thought
    subtime = subage_pos(idx_low:idx_hi)-age_base_pos(i);
    subval_pos(idx_low:idx_hi) = [subtime*5 + val_base_pos(i)] + sd_pos(idx_low:idx_hi);
end

%% Merge simulated positive and negative data into a table
t = table();
t.subid = [subid_neg_tall;subid_pos_tall];
t.age = [subage_neg;subage_pos];
t.val = [subval_neg;subval_pos];
