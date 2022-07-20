function [tout,tdrs] = ILLA(age,value,subid,dt,val0,maxi,skern)
%
% tout = ILLA(age,value,subid,dt,val0,maxi)
% 
% Written By: Tobey J Betthauser, PhD
%             Univsersity of Wisconsin-Madison
%             Department of Medicine
%             Division of Geriatrics
%
% This function is designed to estimate the amyloid vs time curve from
% longitudinal PET imaging data. 

% [msg,id] = lastwarn
warning('off','MATLAB:table:RowsAddedExistingVars');
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
subs = unique(subid);
%% create table object
t = table;
t.age = age;
t.val = value;
t.subid = subid;
t = sortrows(t,{'subid','age'});

%% get within person slopes
t.mrate(:) = nan;
t.max(:) = nan;
t.min(:) = nan;
t.count(:) = nan;
for i = 1:numel(subs)
    idsub = t.subid==subs(i);
    tsub = t(idsub,:);
    if height(tsub)>1
        p = polyfit(tsub.age,tsub.val,1);
        t.mrate(idsub) = p(1);
        t.max(idsub) = max(tsub.val);
        t.min(idsub) = min(tsub.val);
        t.count(idsub) = 1:height(tsub);
    end
    t.nvis(idsub) = height(tsub);
end

%% remove cases with only one scan and simplify table for rate sampling
% tmod = t(~isnan(t.mrate),:);
tmod = t(t.nvis>1,:);
% qval = min(tmod.val):0.01:max(tmod.val);
qval = min(tmod.val):range(tmod.val)/150:max(tmod.val);
tmod = tmod(tmod.count==1,:);
tmod = tmod(tmod.mrate<100,:);

%% Perform descrete rate sampling
tdrs = table();
rate = [];vals = [];
for i = 1:numel(qval)
    ids = tmod.min<qval(i) & tmod.max>qval(i);
%     tdrs.rate(i) = mean(tmod.mrate(ids)); % un-weighted average
    tdrs.val(i) = qval(i);
    tdrs.rate(i) = (tmod.mrate(ids)'*tmod.nvis(ids))/sum(tmod.nvis(ids)); % weighted average
    tdrs.ratestd(i) = std(tmod.mrate(ids));
    tdrs.npos(i) = nnz(tmod.mrate(ids)>0);
    tdrs.tot(i) = nnz(ids);
    
    rate = cat(1,rate,tmod.mrate(ids));
    vals = cat(1,vals,ones(nnz(ids),1)*qval(i));
end
tdrs.ci = 1.96*tdrs.ratestd./sqrt(tdrs.tot); % estimate 95% CI in rate
idkeep = tdrs.tot>=2;
tdrs = tdrs(idkeep,:); % restrict data to observations with more than one observation

% Use robust LOESS to smooth data, unles smoothing kernel is set to null
if skern~=0
    srate = smooth(vals,rate,skern,'rloess');
    [vals,ids] = unique(vals);
    srate = srate(ids);
    tdrs.rate = srate(ismember(vals,tdrs.val));
%     tdrs.rate = smooth(tdrs.val,tdrs.rate,skern,'rloess'); old way of smoothing
end
tdrs.skern(:) = skern;
% Determine if curve is increasing or decreasing with time
med_rate = median(tdrs.rate);

%% Perform iterative model
% set inital conditions
% use iterative process to go forward and backward through descretely
% sampled rate data

% forward in time
qval_cur = mean(tdrs.val);
valf = [];
rf = [];
sdf = [];
nf = [];
nif = 0;
% maxi = 60;
while and(qval_cur<max(tdrs.val),nif<maxi)
    if med_rate<0 && qval_cur<min(tdrs.val)
        break
    end
    
    % below uses standard Euler's Method
    [~,id] = min(abs(tdrs.val - qval_cur)); % lookup slope for current qval
    if tdrs.tot(id)<2
        break
    elseif tdrs.rate(id)<=0 && med_rate>0
        break
    elseif tdrs.rate(id)>=0 && med_rate<0
        break
    end
    valf = cat(1,valf,qval_cur); % store value for current iteration
    nif = nif+1; % iteration counter
    rf = cat(1,rf,tdrs.rate(id)); % store rate for current iteration
    sdf = cat(1,sdf,tdrs.ratestd(id)); % store rate SD for current iteration
    nf = cat(1,nf,tdrs.tot(id)); % store number observations in current iteration
    qval_cur = tdrs.rate(id)*dt + qval_cur; % estimate next querry value using Euler's Method
end
tf = cumsum(dt*ones(nif,1))-dt;

% backward in time
qval_cur = mean(tdrs.val);
valb = [];
rb = [];
sdb = [];
nb = [];
nib = 0;
% maxi = 60;
while and(qval_cur>min(qval),nib<maxi)
    if med_rate<0 && qval_cur>max(tdrs.val)
        break
    end
    [~,id] = min(abs(tdrs.val - qval_cur));
    if tdrs.tot(id)<2
        break
    elseif tdrs.rate(id)<0 && med_rate>0
        break
    elseif tdrs.rate(id)>0 && med_rate<0
        break
    end
    valb = cat(1,valb,qval_cur);
    rb = cat(1,rb,tdrs.rate(id));
    sdb = cat(1,sdb,tdrs.ratestd(id));
    nb = cat(1,nb,tdrs.tot(id));
    qval_cur = tdrs.rate(id)*-dt + qval_cur;
    nib = nib+1;
end
tb = -cumsum(dt*ones(nib,1))+dt;

%% Create output table
tout = table();
tout.val = [flip(valb(2:end));valf];
tout.time = [flip(tb(2:end));tf];
[~,id0] = min(abs(tout.val-val0));
tout.adtime = tout.time-tout.time(id0);
tout.mrate = [flip(rb(2:end));rf];
tout.sdrate = [flip(sdb(2:end));sdf];
tout.nsubs = [flip(nb(2:end));nf];

% Use type I error propagation to estimate model error
tout.sdval = sqrt(((tout.mrate*dt).*sqrt(tout.sdrate.^2./tout.mrate.^2)).^2 + (tout.val*.05).^2); % approximating a 5% error in PET estimates for now
tout.ci95 = 1.96* tout.sdval./sqrt(tout.nsubs);