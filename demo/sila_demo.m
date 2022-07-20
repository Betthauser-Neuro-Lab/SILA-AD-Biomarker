function sila_demo()
% This function was created to show how to train the SILA method and
% estimate time to threshold for individual observations.

% add dependency for simulated data and demo
addpath('demo\')

% simulate data for the example
t = simulate_data();

% train the SILA model
[tsila,tdrs] = SILA(t.age,t.val,t.subid,0.25,21,200);

% Estimate time to threshold for each subject
test = SILA_estimate(tsila,t.age,t.val,t.subid);

% The estimation function has optional input arguments align_years,
% extrap_years, and truncate_aget0. These are input as name-value pairs as
% in the following example:
%   test = SILA_estimate(tsila,t.age,t.val,t.subid,'align_event','all','extrap_years',5,'truncate_aget0','no');

%% These plots show the simulated data and some of the SILA outputs
figure('Units','centimeters','Position',[2,2,12,8])
spaghetti_plot(t.age,t.val,t.subid)
hold on, plot(xlim,21*[1,1],'--k')
title('Simulated Input Data')
xlabel('Age (years)'),ylabel('Value')

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

figure('Units','centimeters','Position',[2,2,12,8])
spaghetti_plot(test.estdtt0,test.val,test.subid)
plot(tsila.adtime,tsila.val,'-k'),hold on
hold on, plot(xlim,21*[1,1],'--k')
title('Data Aligned by Estimated Time to Threshold')
xlabel('Estimated time to threshold (years)'),ylabel('Value')