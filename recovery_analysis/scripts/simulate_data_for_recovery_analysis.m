%Script to simulate data from the three "strong" discrimination models
%Licence: MIT
%Author: Arthur S. Courtin

%% General settings
rng(12345)
addpath(genpath(fullfile(fileparts(cd), 'Palamedes')));

n_datasets = 50;
n_participant = 20;
n_trial = 30;
n_rating_trial = 10;

baseline_temperature = 32;
adapting_temperatures = baseline_temperature + (-2:2);

%% Generate group and participant parameters for the no adaptation model
close all
for n=1:n_datasets
    alpha=-10;
    while any(alpha<-4)
        absolute.dataset(n).mu_alpha=normrnd(0,1);
        absolute.dataset(n).mu_log_beta=normrnd(0,1);
        absolute.dataset(n).mu_hlogit_lambda=normrnd(-4,1);

        absolute.dataset(n).mu_intercept=normrnd(-2,1);
        absolute.dataset(n).mu_log_slope=normrnd(-2,1);
        absolute.dataset(n).mu_log_lb=normrnd(2,.5);
        absolute.dataset(n).mu_log_ub=normrnd(2,.5);
        absolute.dataset(n).mu_log_eta=normrnd(3,1);

        absolute.dataset(n).tau_alpha=abs(normrnd(0,1));
        absolute.dataset(n).tau_log_beta=abs(normrnd(0,1));
        absolute.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));

        absolute.dataset(n).tau_intercept=abs(normrnd(-2,1));
        absolute.dataset(n).tau_log_slope=abs(normrnd(-2,1));
        absolute.dataset(n).tau_log_lb=abs(normrnd(2,.5));
        absolute.dataset(n).tau_log_ub=abs(normrnd(2,.5));
        absolute.dataset(n).tau_log_eta=abs(normrnd(3,1));

        for p=1:n_participant
            alpha(p)=normrnd(absolute.dataset(n).mu_alpha,absolute.dataset(n).tau_alpha);
            absolute.dataset(n).participant(p).alpha=alpha(p);
            absolute.dataset(n).participant(p).beta=exp(normrnd(absolute.dataset(n).mu_log_beta,absolute.dataset(n).tau_log_beta));
            absolute.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(absolute.dataset(n).mu_hlogit_lambda,absolute.dataset(n).tau_hlogit_lambda)));

            absolute.dataset(n).participant(p).intercept=normrnd(absolute.dataset(n).mu_intercept,absolute.dataset(n).tau_intercept);
            absolute.dataset(n).participant(p).slope=exp(normrnd(absolute.dataset(n).mu_log_slope,absolute.dataset(n).tau_log_slope));
            absolute.dataset(n).participant(p).lb=exp(normrnd(absolute.dataset(n).mu_log_lb,absolute.dataset(n).tau_log_lb));
            absolute.dataset(n).participant(p).ub=exp(normrnd(absolute.dataset(n).mu_log_ub,absolute.dataset(n).tau_log_ub));
            absolute.dataset(n).participant(p).eta=exp(normrnd(absolute.dataset(n).mu_log_eta,absolute.dataset(n).tau_log_eta));
        end
    end
%     x=30:.1:40;
%     figure()
%     hold on
%     for p=1:n_participant
%         x_c=x-34-absolute.dataset(n).participant(p).alpha;
%         weight = inv_logit(x_c.*100);
%         theta = (1-weight) * 0.5 + (weight-0) .* normcdf(absolute.dataset(n).participant(p).beta .* x_c);
%         plot(x,theta)
%     end
%     hold off
%     xlim([30 40])
%     ylim([.5 1])
end
%% Generate group and participant parameters for the full adaptation model
close all
for n=1:n_datasets
    relative.dataset(n).mu_log_alpha=normrnd(-2,1);
    relative.dataset(n).mu_log_beta=normrnd(0,1);
    relative.dataset(n).mu_hlogit_lambda=normrnd(-4,1);
    
    relative.dataset(n).mu_intercept=normrnd(-2,1);
    relative.dataset(n).mu_log_slope=normrnd(-2,1);
    relative.dataset(n).mu_log_lb=normrnd(2,.5);
    relative.dataset(n).mu_log_ub=normrnd(2,.5);
    relative.dataset(n).mu_log_eta=normrnd(3,1);
    
    relative.dataset(n).tau_log_alpha=abs(normrnd(0,1));
    relative.dataset(n).tau_log_beta=abs(normrnd(0,1));
    relative.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));
    
    relative.dataset(n).tau_intercept=abs(normrnd(-2,1));
    relative.dataset(n).tau_log_slope=abs(normrnd(-2,1));
    relative.dataset(n).tau_log_lb=abs(normrnd(2,.5));
    relative.dataset(n).tau_log_ub=abs(normrnd(2,.5));
    relative.dataset(n).tau_log_eta=abs(normrnd(3,1));
    
    for p=1:n_participant
        relative.dataset(n).participant(p).alpha=exp(normrnd(relative.dataset(n).mu_log_alpha,relative.dataset(n).tau_log_alpha));
        relative.dataset(n).participant(p).beta=exp(normrnd(relative.dataset(n).mu_log_beta,relative.dataset(n).tau_log_beta));
        relative.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(relative.dataset(n).mu_hlogit_lambda,relative.dataset(n).tau_hlogit_lambda)));
        
        relative.dataset(n).participant(p).intercept=normrnd(relative.dataset(n).mu_intercept,relative.dataset(n).tau_intercept);
        relative.dataset(n).participant(p).slope=exp(normrnd(relative.dataset(n).mu_log_slope,relative.dataset(n).tau_log_slope));
        relative.dataset(n).participant(p).lb=exp(normrnd(relative.dataset(n).mu_log_lb,relative.dataset(n).tau_log_lb));
        relative.dataset(n).participant(p).ub=exp(normrnd(relative.dataset(n).mu_log_ub,relative.dataset(n).tau_log_ub));
        relative.dataset(n).participant(p).eta=exp(normrnd(relative.dataset(n).mu_log_eta,relative.dataset(n).tau_log_eta));
    end
    
%     x=0:.1:10;
%     figure()
%     hold on
%     for p=1:n_participant
%         x_c=x-relative.dataset(n).participant(p).alpha;
%         weight = inv_logit(x_c.*100);
%         theta = (1-weight) * 0.5 + (weight-0) .* normcdf(relative.dataset(n).participant(p).beta .* x_c);
%         plot(x,theta)
%     end
%     hold off
%     xlim([0 10])
%     ylim([.5 1])
end
%% Generate group and participant parameters for the partial adaptation model
close all
for n=1:n_datasets
    sat=1000;
    while any((mean(sat)>45)+(sum(sat>50)))
        mixed.dataset(n).mu_log_alpha=normrnd(-2,1);
        mixed.dataset(n).mu_log_sigma=normrnd(1,1);
        mixed.dataset(n).mu_hlogit_lambda=normrnd(-4,1);
        
        mixed.dataset(n).mu_intercept=normrnd(-2,1);
        mixed.dataset(n).mu_log_slope=normrnd(-2,1);
        mixed.dataset(n).mu_log_lb=normrnd(2,.5);
        mixed.dataset(n).mu_log_ub=normrnd(2,.5);
        mixed.dataset(n).mu_log_eta=normrnd(3,1);

        mixed.dataset(n).tau_log_alpha=abs(normrnd(0,1));
        mixed.dataset(n).tau_log_sigma=abs(normrnd(0,1));
        mixed.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));
        
        mixed.dataset(n).tau_intercept=abs(normrnd(-2,1));
        mixed.dataset(n).tau_log_slope=abs(normrnd(-2,1));
        mixed.dataset(n).tau_log_lb=abs(normrnd(2,.5));
        mixed.dataset(n).tau_log_ub=abs(normrnd(2,.5));
        mixed.dataset(n).tau_log_eta=abs(normrnd(3,1));

        for p=1:n_participant
            alpha=exp(normrnd(mixed.dataset(n).mu_log_alpha,mixed.dataset(n).tau_log_alpha));
            sigma=exp(normrnd(mixed.dataset(n).mu_log_sigma,mixed.dataset(n).tau_log_sigma));
            mixed.dataset(n).participant(p).alpha=alpha;
            mixed.dataset(n).participant(p).sigma=sigma;
            mixed.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(mixed.dataset(n).mu_hlogit_lambda,mixed.dataset(n).tau_hlogit_lambda)));
        
            mixed.dataset(n).participant(p).intercept=normrnd(mixed.dataset(n).mu_intercept,mixed.dataset(n).tau_intercept);
            mixed.dataset(n).participant(p).slope=exp(normrnd(mixed.dataset(n).mu_log_slope,mixed.dataset(n).tau_log_slope));
            mixed.dataset(n).participant(p).lb=exp(normrnd(mixed.dataset(n).mu_log_lb,mixed.dataset(n).tau_log_lb));
            mixed.dataset(n).participant(p).ub=exp(normrnd(mixed.dataset(n).mu_log_ub,mixed.dataset(n).tau_log_ub));
            mixed.dataset(n).participant(p).eta=exp(normrnd(mixed.dataset(n).mu_log_eta,mixed.dataset(n).tau_log_eta));
        end
        for p=1:n_participant
            sat(p)=34+mixed.dataset(n).participant(p).sigma+mixed.dataset(n).participant(p).alpha;
        end
    end
%     x=30:.2:50;
%     figure()
%     hold on
%     for p=1:n_participant
%         for at=adapting_temperatures
%             x_r=x-at;
%             x_c=x_r-mixed.dataset(n).participant(p).alpha;
%             beta=3/(2+mixed.dataset(n).participant(p).sigma+mixed.dataset(n).participant(p).alpha-(at-32));
%             weight = inv_logit(x_c.*100);
%             theta = (1-weight) * 0.5 + (weight-0) .* normcdf(beta .* x_c);
%             plot(x,theta)
%         end
%     end
%     hold off
%     ylim([.5 1])
%     s(n,1:n_participant)=sat;
end
% plot(s,'o')
% hold on
% plot([1 50],45*ones(1,2))
% hold off
%% Initialize  generic PM (matching warm discrimination settings)
close all
PM = PAL_AMPM_setupPM( ...
    'priorAlphaRange',0.2:0.2:10,...
    'priorBetaRange',0:0.05:1.5,...
    'priorLambdaRange',0.01:0.02:0.19,...
    'priorGammaRange',0.5,...
    'stimRange',0.2:0.2:10,...
    'PF',@PAL_Quick...
    );
prior = ...
    PAL_pdfNormal(PM.priorAlphas,2,2).*...
    ones(size(PM.priorBetas)).*...
    ones(size(PM.priorGammas)).*...
    PAL_pdfNormal(log(2*PM.priorLambdas./(1-2*PM.priorLambdas)),-4,1);

prior=prior./sum(sum(sum(sum(prior))));
    
PM = PAL_AMPM_setupPM(PM,'prior',prior);
%% Simulate trial data for the absolute model
for n=1:n_datasets
    fprintf('Simulating absolute coding dataset %i out of %i \n',n,n_datasets)
    for p=1:n_participant
        fprintf('Participant %i out of %i \n',p,n_participant)
        alpha=absolute.dataset(n).participant(p).alpha;
        beta=absolute.dataset(n).participant(p).beta;
        lambda=absolute.dataset(n).participant(p).lambda;        
        
        intercept=absolute.dataset(n).participant(p).intercept;
        slope=absolute.dataset(n).participant(p).slope;
        lb=absolute.dataset(n).participant(p).lb;
        ub=absolute.dataset(n).participant(p).ub;
        eta=absolute.dataset(n).participant(p).eta;

        for at=adapting_temperatures
            adx=at-32+3;
            PMl=PM;
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
                centered_stimulus = absolute_target - 34 - alpha;
                weight = inv_logit(centered_stimulus.*100);
                theta = (1-weight) * 0.5 + (weight-lambda) .* normcdf(beta .* centered_stimulus);
                
                choice_accuracy=binornd(1,theta);
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                absolute.dataset(n).participant(p).condition(adx).baseline_temperature=baseline_temperature;
                absolute.dataset(n).participant(p).condition(adx).absolute_adapting_temperature=at;
                
                absolute.dataset(n).participant(p).condition(adx).trial(t)=t;
                absolute.dataset(n).participant(p).condition(adx).absolute_target_temperature(t)=absolute_target;
                absolute.dataset(n).participant(p).condition(adx).choice_accuracy(t)=choice_accuracy;
            end
            
            [target_temperatures] = determine_target_detection_at(PMl,at);
            for t=1:n_rating_trial
                if mod(t,2)
                    tt=target_temperatures(1);
                else
                    tt=target_temperatures(2);
                end
                lr=intercept+slope*(tt-42);

                r0=binornd(1,1-inv_logit(lr+lb));
                r1=binornd(1,inv_logit(lr-ub));
                mr=inv_logit(lr);
                r=betarnd(mr*eta,(1-mr)*eta);
                
                if r0
                    rating=0;
                elseif r1
                    rating=1;
                else
                    rating=r;
                end

                absolute.dataset(n).participant(p).condition(adx).rating_trial(t)=t;
                absolute.dataset(n).participant(p).condition(adx).rating_absolute_target_temperature(t)=tt;
                absolute.dataset(n).participant(p).condition(adx).rating(t)=rating;
            end
        end
    end
end
%% Reformat and save the absolute coding discrimination simulations as csv
rows={};
row_idx=1;
for n=1:n_datasets
    
    ma=absolute.dataset(n).mu_alpha;
    mlb=absolute.dataset(n).mu_log_beta;
    mll=absolute.dataset(n).mu_hlogit_lambda;
    
    ta=absolute.dataset(n).tau_alpha;
    tlb=absolute.dataset(n).tau_log_beta;
    tll=absolute.dataset(n).tau_hlogit_lambda; 
        
    for p=1:n_participant
        alpha=absolute.dataset(n).participant(p).alpha;
        beta=absolute.dataset(n).participant(p).beta;
        lambda=absolute.dataset(n).participant(p).lambda;        

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_trial
                tt=absolute.dataset(n).participant(p).condition(adx).absolute_target_temperature(t);
                ca=absolute.dataset(n).participant(p).condition(adx).choice_accuracy(t);
                    
                rows(row_idx,:)={...
                    n,'a',...
                    ma,mlb,mll,...
                    ta,tlb,tll,...
                    p,...
                    alpha,beta,lambda,...
                    baseline_temperature,at,...
                    t,tt,ca...
                    };
                row_idx=row_idx+1;
            end
        end
    end
end
absolute_table = cell2table(rows, ...
    'VariableNames', { ...
        'dataset', ...
        'model', ...
        'mu_alpha',...
        'mu_log_beta',...
        'mu_hlogit_lambda',...
        'tau_alpha',...
        'tau_log_beta',...
        'tau_hlogit_lambda',...
        'participant', ...
        'alpha', ...
        'beta', ...
        'lambda', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'choice_accuracy' ...
    });
outputPath = fullfile(fileparts(cd), 'simulated_data', 'absolute_model_discrimination_data.csv');
writetable(absolute_table, outputPath);
%% Reformat and save the absolute coding rating simulations as csv
rows={};
row_idx=1;
for n=1:n_datasets
    
    mi=absolute.dataset(n).mu_intercept;
    mls=absolute.dataset(n).mu_log_slope;
    mllb=absolute.dataset(n).mu_log_lb;
    mlub=absolute.dataset(n).mu_log_ub;
    mle=absolute.dataset(n).mu_log_eta;

    ti=absolute.dataset(n).tau_intercept;
    tls=absolute.dataset(n).tau_log_slope;
    tllb=absolute.dataset(n).tau_log_lb;
    tlub=absolute.dataset(n).tau_log_ub;
    tle=absolute.dataset(n).tau_log_eta; 
        
    for p=1:n_participant
        intercept=absolute.dataset(n).participant(p).intercept;
        slope=absolute.dataset(n).participant(p).slope;
        lb=absolute.dataset(n).participant(p).lb;        
        ub=absolute.dataset(n).participant(p).ub;        
        eta=absolute.dataset(n).participant(p).eta;        

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_rating_trial
                tt=absolute.dataset(n).participant(p).condition(adx).rating_absolute_target_temperature(t);
                r=absolute.dataset(n).participant(p).condition(adx).rating(t);
                    
                rows(row_idx,:)={...
                    n,'a',...
                    mi,mls,mllb,mlub,mle,...
                    ti,tls,tllb,tlub,tle,...
                    p,...
                    intercept,slope,lb,ub,eta,...
                    baseline_temperature,at,...
                    t,tt,r...
                    };
                row_idx=row_idx+1;
            end
        end
    end
end
absolute_table = cell2table(rows, ...
    'VariableNames', { ...
        'dataset', ...
        'model', ...
        'mu_intercept',...
        'mu_log_slope',...
        'mu_log_lb',...
        'mu_log_ub',...
        'mu_log_eta',...
        'tau_intercept',...
        'tau_log_slope',...
        'tau_log_lb',...
        'tau_log_ub',...
        'tau_log_eta',...
        'participant', ...
        'intercept', ...
        'slope', ...
        'lb', ...
        'ub', ...
        'eta', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'rating' ...
    });
outputPath = fullfile(fileparts(cd), 'simulated_data', 'absolute_model_rating_data.csv');
writetable(absolute_table, outputPath);

%% Simulate trial data for the relative model
for n=1:n_datasets
    fprintf('Simulating relative coding dataset %i out of %i \n',n,n_datasets)
    for p=1:n_participant
        fprintf('Participant %i out of %i \n',p,n_participant)
        alpha=relative.dataset(n).participant(p).alpha;
        beta=relative.dataset(n).participant(p).beta;
        lambda=relative.dataset(n).participant(p).lambda;        

        intercept=relative.dataset(n).participant(p).intercept;
        slope=relative.dataset(n).participant(p).slope;
        lb=relative.dataset(n).participant(p).lb;
        ub=relative.dataset(n).participant(p).ub;
        eta=relative.dataset(n).participant(p).eta;
        
        for at=adapting_temperatures
            adx=at-32+3;
            PMl=PM;
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
                
                centered_stimulus = relative_target - alpha;
                weight = inv_logit(centered_stimulus.*100);
                theta = (1-weight) * 0.5 + (weight-lambda) .* normcdf(beta .* centered_stimulus);
                
                choice_accuracy=binornd(1,theta);
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                relative.dataset(n).participant(p).condition(adx).baseline_temperature=baseline_temperature;
                relative.dataset(n).participant(p).condition(adx).absolute_adapting_temperature=at;
                
                relative.dataset(n).participant(p).condition(adx).trial(t)=t;
                relative.dataset(n).participant(p).condition(adx).absolute_target_temperature(t)=absolute_target;
                relative.dataset(n).participant(p).condition(adx).choice_accuracy(t)=choice_accuracy;
            end
            
            [target_temperatures] = determine_target_detection_at(PMl,at);
            for t=1:n_rating_trial
                if mod(t,2)
                    tt=target_temperatures(1);
                else
                    tt=target_temperatures(2);
                end
                lr=intercept+slope*(tt-at-8);

                r0=binornd(1,1-inv_logit(lr+lb));
                r1=binornd(1,inv_logit(lr-ub));
                mr=inv_logit(lr);
                r=betarnd(mr*eta,(1-mr)*eta);
                
                if r0
                    rating=0;
                elseif r1
                    rating=1;
                else
                    rating=r;
                end

                relative.dataset(n).participant(p).condition(adx).rating_trial(t)=t;
                relative.dataset(n).participant(p).condition(adx).rating_absolute_target_temperature(t)=tt;
                relative.dataset(n).participant(p).condition(adx).rating(t)=rating;
            end
        end
    end
end
%% Reformat and save the relative coding simulations as csv
rows={};
row_idx=1;
for n=1:n_datasets
    
    mla=relative.dataset(n).mu_log_alpha;
    mlb=relative.dataset(n).mu_log_beta;
    mll=relative.dataset(n).mu_hlogit_lambda;
    
    tla=relative.dataset(n).tau_log_alpha;
    tlb=relative.dataset(n).tau_log_beta;
    tll=relative.dataset(n).tau_hlogit_lambda; 
        
    for p=1:n_participant
        alpha=relative.dataset(n).participant(p).alpha;
        beta=relative.dataset(n).participant(p).beta;
        lambda=relative.dataset(n).participant(p).lambda;        

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_trial
                tt=relative.dataset(n).participant(p).condition(adx).absolute_target_temperature(t);
                ca=relative.dataset(n).participant(p).condition(adx).choice_accuracy(t);
                    
                rows(row_idx,:)={...
                    n,'r',...
                    mla,mlb,mll,...
                    tla,tlb,tll,...
                    p,...
                    alpha,beta,lambda,...
                    baseline_temperature,at,...
                    t,tt,ca...
                    };
                row_idx=row_idx+1;
            end
        end
    end
end
relative_table = cell2table(rows, ...
    'VariableNames', { ...
        'dataset', ...
        'model', ...
        'mu_log_alpha',...
        'mu_log_beta',...
        'mu_hlogit_lambda',...
        'tau_log_alpha',...
        'tau_log_beta',...
        'tau_hlogit_lambda',...
        'participant', ...
        'alpha', ...
        'beta', ...
        'lambda', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'choice_accuracy' ...
    });
outputPath = fullfile(fileparts(cd), 'simulated_data', 'relative_model_discrimination_data.csv');
writetable(relative_table, outputPath);

%% Reformat and save the relative coding rating simulations as csv
rows={};
row_idx=1;
for n=1:n_datasets
    
    mi=relative.dataset(n).mu_intercept;
    mls=relative.dataset(n).mu_log_slope;
    mllb=relative.dataset(n).mu_log_lb;
    mlub=relative.dataset(n).mu_log_ub;
    mle=relative.dataset(n).mu_log_eta;

    ti=relative.dataset(n).tau_intercept;
    tls=relative.dataset(n).tau_log_slope;
    tllb=relative.dataset(n).tau_log_lb;
    tlub=relative.dataset(n).tau_log_ub;
    tle=relative.dataset(n).tau_log_eta; 
        
    for p=1:n_participant
        intercept=relative.dataset(n).participant(p).intercept;
        slope=relative.dataset(n).participant(p).slope;
        lb=relative.dataset(n).participant(p).lb;        
        ub=relative.dataset(n).participant(p).ub;        
        eta=relative.dataset(n).participant(p).eta;        

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_rating_trial
                tt=relative.dataset(n).participant(p).condition(adx).rating_absolute_target_temperature(t);
                r=relative.dataset(n).participant(p).condition(adx).rating(t);
                    
                rows(row_idx,:)={...
                    n,'r',...
                    mi,mls,mllb,mlub,mle,...
                    ti,tls,tllb,tlub,tle,...
                    p,...
                    intercept,slope,lb,ub,eta,...
                    baseline_temperature,at,...
                    t,tt,r...
                    };
                row_idx=row_idx+1;
            end
        end
    end
end
relative_table = cell2table(rows, ...
    'VariableNames', { ...
        'dataset', ...
        'model', ...
        'mu_intercept',...
        'mu_log_slope',...
        'mu_log_lb',...
        'mu_log_ub',...
        'mu_log_eta',...
        'tau_intercept',...
        'tau_log_slope',...
        'tau_log_lb',...
        'tau_log_ub',...
        'tau_log_eta',...
        'participant', ...
        'intercept', ...
        'slope', ...
        'lb', ...
        'ub', ...
        'eta', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'rating' ...
    });
outputPath = fullfile(fileparts(cd), 'simulated_data', 'relative_model_rating_data.csv');
writetable(relative_table, outputPath);

%% Simulate trial data for the mixed model
for n=1:n_datasets
    fprintf('Simulating mixed coding dataset %i out of %i \n',n,n_datasets)
    for p=1:n_participant
        fprintf('Participant %i out of %i \n',p,n_participant)
        alpha=mixed.dataset(n).participant(p).alpha;
        sigma=mixed.dataset(n).participant(p).sigma;
        lambda=mixed.dataset(n).participant(p).lambda;        

        intercept=mixed.dataset(n).participant(p).intercept;
        slope=mixed.dataset(n).participant(p).slope;
        lb=mixed.dataset(n).participant(p).lb;
        ub=mixed.dataset(n).participant(p).ub;
        eta=mixed.dataset(n).participant(p).eta;
        
        for at=adapting_temperatures
            adx=at-32+3;
            PMl=PM;
            beta = 3/(34+sigma+alpha-at);
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
           
                centered_stimulus = relative_target - alpha;
                weight = inv_logit(centered_stimulus.*100);
                theta = (1-weight) * 0.5 + (weight-lambda) .* normcdf(beta .* centered_stimulus);
                
                choice_accuracy=binornd(1,theta);
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                mixed.dataset(n).participant(p).condition(adx).baseline_temperature=baseline_temperature;
                mixed.dataset(n).participant(p).condition(adx).absolute_adapting_temperature=at;
                
                mixed.dataset(n).participant(p).condition(adx).trial(t)=t;
                mixed.dataset(n).participant(p).condition(adx).absolute_target_temperature(t)=absolute_target;
                mixed.dataset(n).participant(p).condition(adx).choice_accuracy(t)=choice_accuracy;
            end
            
            [target_temperatures] = determine_target_detection_at(PMl,at);
            for t=1:n_rating_trial
                if mod(t,2)
                    tt=target_temperatures(1);
                else
                    tt=target_temperatures(2);
                end
                lr=intercept+slope*(tt-at-8)*beta;

                r0=binornd(1,1-inv_logit(lr+lb));
                r1=binornd(1,inv_logit(lr-ub));
                mr=inv_logit(lr);
                r=betarnd(mr*eta,(1-mr)*eta);
                
                if r0
                    rating=0;
                elseif r1
                    rating=1;
                else
                    rating=r;
                end

                mixed.dataset(n).participant(p).condition(adx).rating_trial(t)=t;
                mixed.dataset(n).participant(p).condition(adx).rating_absolute_target_temperature(t)=tt;
                mixed.dataset(n).participant(p).condition(adx).rating(t)=rating;
            end
        end
    end
end

%% Reformat and save the mixed coding simulations as csv
rows={};
row_idx=1;
for n=1:n_datasets
    
    mla=mixed.dataset(n).mu_log_alpha;
    mls=mixed.dataset(n).mu_log_sigma;
    mll=mixed.dataset(n).mu_hlogit_lambda;
    
    tla=mixed.dataset(n).tau_log_alpha;
    tls=mixed.dataset(n).tau_log_sigma;
    tll=mixed.dataset(n).tau_hlogit_lambda; 
        
    for p=1:n_participant
        alpha=mixed.dataset(n).participant(p).alpha;
        sigma=mixed.dataset(n).participant(p).sigma;
        lambda=mixed.dataset(n).participant(p).lambda;        

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_trial
                tt=mixed.dataset(n).participant(p).condition(adx).absolute_target_temperature(t);
                ca=mixed.dataset(n).participant(p).condition(adx).choice_accuracy(t);
                    
                rows(row_idx,:)={...
                    n,'m',...
                    mla,mls,mll,...
                    tla,tls,tll,...
                    p,...
                    alpha,sigma,lambda,...
                    baseline_temperature,at,...
                    t,tt,ca...
                    };
                row_idx=row_idx+1;
            end
        end
    end
end
mixed_table = cell2table(rows, ...
    'VariableNames', { ...
        'dataset', ...
        'model', ...
        'mu_log_alpha',...
        'mu_log_sigma',...
        'mu_hlogit_lambda',...
        'tau_log_alpha',...
        'tau_log_sigma',...
        'tau_hlogit_lambda',...
        'participant', ...
        'alpha', ...
        'sigma', ...
        'lambda', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'choice_accuracy' ...
    });
outputPath = fullfile(fileparts(cd), 'simulated_data', 'mixed_model_discrimination_data.csv');
writetable(mixed_table, outputPath);

%% Reformat and save the mixed coding rating simulations as csv
rows={};
row_idx=1;
for n=1:n_datasets
    
    mi=mixed.dataset(n).mu_intercept;
    mls=mixed.dataset(n).mu_log_slope;
    mllb=mixed.dataset(n).mu_log_lb;
    mlub=mixed.dataset(n).mu_log_ub;
    mle=mixed.dataset(n).mu_log_eta;

    ti=mixed.dataset(n).tau_intercept;
    tls=mixed.dataset(n).tau_log_slope;
    tllb=mixed.dataset(n).tau_log_lb;
    tlub=mixed.dataset(n).tau_log_ub;
    tle=mixed.dataset(n).tau_log_eta; 
        
    for p=1:n_participant
        intercept=mixed.dataset(n).participant(p).intercept;
        slope=mixed.dataset(n).participant(p).slope;
        lb=mixed.dataset(n).participant(p).lb;        
        ub=mixed.dataset(n).participant(p).ub;        
        eta=mixed.dataset(n).participant(p).eta;        

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_rating_trial
                tt=mixed.dataset(n).participant(p).condition(adx).rating_absolute_target_temperature(t);
                r=mixed.dataset(n).participant(p).condition(adx).rating(t);
                    
                rows(row_idx,:)={...
                    n,'m',...
                    mi,mls,mllb,mlub,mle,...
                    ti,tls,tllb,tlub,tle,...
                    p,...
                    intercept,slope,lb,ub,eta,...
                    baseline_temperature,at,...
                    t,tt,r...
                    };
                row_idx=row_idx+1;
            end
        end
    end
end
mixed_table = cell2table(rows, ...
    'VariableNames', { ...
        'dataset', ...
        'model', ...
        'mu_intercept',...
        'mu_log_slope',...
        'mu_log_lb',...
        'mu_log_ub',...
        'mu_log_eta',...
        'tau_intercept',...
        'tau_log_slope',...
        'tau_log_lb',...
        'tau_log_ub',...
        'tau_log_eta',...
        'participant', ...
        'intercept', ...
        'slope', ...
        'lb', ...
        'ub', ...
        'eta', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'rating' ...
    });
outputPath = fullfile(fileparts(cd), 'simulated_data', 'mixed_model_rating_data.csv');
writetable(mixed_table, outputPath);


%%
function y=inv_logit(x) 
    y=1./(1+exp(-x));
end

function [target_temperatures] = determine_target_detection_at(PM,at)
    x = PM.stimRange;
    
    for idx=1:length(x)
        y(idx)=sum(PM.response(PM.x(1:end-1)==x(idx)));
        n(idx)=sum(PM.x(1:end-1)==x(idx));
    end
    
    grid.alpha=PM.priorAlphaRange;
    grid.beta=10.^PM.priorBetaRange;
    grid.gamma=PM.priorGammaRange;
    grid.lambda=PM.priorLambdaRange;
    
    [paramsValues, LL, scenario, output] = PAL_PFML_Fit(x, y, n, grid, [1 1 0 1], @PAL_Quick);
    paramsValues(4)=0;
    
    if scenario==1
    else
        paramsValues(1)=PM.threshold(end);
        paramsValues(2)=10^(PM.slope(end));
    end
    
    target_temperatures = [at+PAL_Quick(paramsValues, 1-10^-3, 'Inverse') at+8];
  
end

