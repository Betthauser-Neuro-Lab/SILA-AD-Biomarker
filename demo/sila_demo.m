function sila_demo()
% This function was created to show how to train the SILA method and
% estimate time to threshold for individual observations.

%% add filepath for SILA functions to Matlab search path
path_demo = fileparts(mfilename('fullpath'));
addpath(fullfile(path_demo,'..'))

%% simulate data for the example
disp('Simulating Data')
t = simulate_data();

%% train the SILA model
disp('Training the SILA model using SILA.m')
[tsila,tdrs] = SILA(t.age,t.val,t.subid,0.25,21,200);

% in this demo, we are inputting tall-format age, val, and subid and
% specifying to use 0.25 year intervals for numeric integration with a
% threshold value of 21 and a maximum number of iterations from the median
% of 200. With these inputs, the maximum amount of time that can be modeled
% is 200 iterations x 0.25 years/iteration = 50 years in each direction.

%% Estimate time to threshold and age at threshold for each subject
disp('Generating subject-level estimates with SILA_estimate.m')
test = SILA_estimate(tsila,t.age,t.val,t.subid);

% The estimation function has optional input arguments align_years,
% extrap_years, and truncate_aget0. These are input as name-value pairs as
% in the following example:
%   test = SILA_estimate(tsila,t.age,t.val,t.subid,'align_event','all','extrap_years',5,'truncate_aget0','no');

%% These plots show the simulated data and some of the SILA outputs
disp('Generating plots of the data')

% spaghetti plot of value vs. age for simulated data
figure('Units','centimeters','Position',[2,2,12,8])
spaghetti_plot(t.age,t.val,t.subid)
hold on, plot(xlim,21*[1,1],'--k')
title('Simulated Input Data')
xlabel('Age (years)'),ylabel('Value')

% plots showing the output from descrete rate sampling (i.e., rate vs. value) 
% and modeled value vs. time data.
figure('Units','centimeters','Position',[2,2,12,12])
subplot(2,1,1)
plot(tdrs.val,tdrs.rate,'-'),hold on
plot(tdrs.val,tdrs.rate + tdrs.ci,'--r')
plot(tdrs.val,tdrs.rate - tdrs.ci,'--r')
title('Discrete Rate Sampling Curve')
xlabel('Value'),ylabel('\DeltaValue per Year')

subplot(2,1,2)
plot(tsila.adtime,tsila.val,'-'),hold on
plot(xlim,21*[1,1],'--k')
title('SILA Modeled{\it Value vs. Time} Curve')
xlabel('Time from Threshold'),ylabel('Value')
legend({'Modeled curve','threshold'},'Location','northwest')

% value vs. time for all subjects
figure('Units','centimeters','Position',[2,2,12,8])
spaghetti_plot(test.estdtt0,test.val,test.subid)
plot(tsila.adtime,tsila.val,'-k'),hold on
hold on, plot(xlim,21*[1,1],'--k')
title('Data Aligned by Estimated Time to Threshold')
xlabel('Estimated time to threshold (years)'),ylabel('Value')

% value vs. time for an indivdual case
sub = find(test.estdtt0>1 & test.estdtt0<10,1);
ids = test.subid==test.subid(sub);

figure('Units','centimeters','Position',[2,2,9,12])
subplot(2,1,1)
spaghetti_plot(test.age(ids),test.val(ids),test.subid(ids))
hold on, plot([min(t.age),max(t.age)],21*[1,1],'--k')
title('Observations by Age')
xlabel('Age (years)'),ylabel('Value')
legend({'Individual Case Observations'},'Location','northwest')
ylim([min(tsila.val),max(tsila.val)])
xlim([min(t.age),max(t.age)])

subplot(2,1,2)
spaghetti_plot(test.estdtt0(ids),test.val(ids),test.subid(ids))
plot(tsila.adtime,tsila.val,'-k'),hold on
hold on, plot(xlim,21*[1,1],'--k')
xlim([min(tsila.adtime),max(tsila.adtime)])
title('Observations by Estimated Time to Threshold')
xlabel('Estimated time to threshold (years)'),ylabel('Value')
legend({'Individual Case Observations','SILA Modeled Values'},'Location','northwest')