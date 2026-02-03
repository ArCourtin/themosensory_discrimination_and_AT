%Script to simulate data from the three "strong" discrimination models
%Licence: MIT
%Author: Arthur S. Courtin

%% General settings
rng(12345)
addpath(genpath([cd '\Palameses\']));

n_datasets = 50;
n_participant = 20;
n_trial = 30;

baseline_temperature = 32;
adapting_temperatures = baseline_temperature + (-2:2);

%% Generate group and participant parameters for the no adaptation model
close all
for n=1:n_datasets
    absolute.dataset(n).mu_log_alpha=normrnd(-2,1);
    absolute.dataset(n).mu_log_beta=normrnd(0,1);
    absolute.dataset(n).mu_hlogit_lambda=normrnd(-4,1);
    
    absolute.dataset(n).tau_log_alpha=abs(normrnd(0,1));
    absolute.dataset(n).tau_log_beta=abs(normrnd(0,1));
    absolute.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));
    
    for p=1:n_participant
        absolute.dataset(n).participant(p).alpha=exp(normrnd(absolute.dataset(n).mu_log_alpha,absolute.dataset(n).tau_log_alpha));
        absolute.dataset(n).participant(p).beta=exp(normrnd(absolute.dataset(n).mu_log_beta,absolute.dataset(n).tau_log_beta));
        absolute.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(absolute.dataset(n).mu_hlogit_lambda,absolute.dataset(n).tau_hlogit_lambda)));
    end
    
%     x=30:.1:40;
%     figure()
%     hold on
%     for p=1:n_participant
%         plot(x,normcdf(x,absolute.dataset(n).participant(p).alpha+baseline_temperature+2,1/absolute.dataset(n).participant(p).beta))
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
    
    relative.dataset(n).tau_log_alpha=abs(normrnd(0,1));
    relative.dataset(n).tau_log_beta=abs(normrnd(0,1));
    relative.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));
    
    for p=1:n_participant
        relative.dataset(n).participant(p).alpha=exp(normrnd(relative.dataset(n).mu_log_alpha,relative.dataset(n).tau_log_alpha));
        relative.dataset(n).participant(p).beta=exp(normrnd(relative.dataset(n).mu_log_beta,relative.dataset(n).tau_log_beta));
        relative.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(relative.dataset(n).mu_hlogit_lambda,relative.dataset(n).tau_hlogit_lambda)));
    end
    
%     x=0:.1:10;
%     figure()
%     hold on
%     for p=1:n_participant
%         plot(x,normcdf(x,relative.dataset(n).participant(p).alpha,1/relative.dataset(n).participant(p).beta))
%     end
%     hold off
%     xlim([0 10])
%     ylim([.5 1])
end
%% Generate group and participant parameters for the partial adaptation model
close all
for n=1:n_datasets
    mixed.dataset(n).mu_log_alpha=normrnd(-2,1);
    mixed.dataset(n).mu_log_sigma=normrnd(0,1);
    mixed.dataset(n).mu_hlogit_lambda=normrnd(-4,1);
    
    mixed.dataset(n).tau_log_alpha=abs(normrnd(0,1));
    mixed.dataset(n).tau_log_sigma=abs(normrnd(0,1));
    mixed.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));
    
    for p=1:n_participant
        alpha=exp(normrnd(mixed.dataset(n).mu_log_alpha,mixed.dataset(n).tau_log_alpha));
        sigma=exp(normrnd(mixed.dataset(n).mu_log_sigma,mixed.dataset(n).tau_log_sigma));
        mixed.dataset(n).participant(p).alpha=alpha;
        mixed.dataset(n).participant(p).sigma=sigma;
        mixed.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(mixed.dataset(n).mu_hlogit_lambda,mixed.dataset(n).tau_hlogit_lambda)));
    end
%     
%     x=-2:.1:10;
%     figure()
%     hold on
%     for p=1:n_participant
%         for at=adapting_temperatures
%             beta=3/(mixed.dataset(n).participant(p).sigma+2-(at-32));
%             plot(x,normcdf(beta*(x-mixed.dataset(n).participant(p).alpha-(at-32))))
%         end
%     end
%     hold off
%     ylim([.5 1])
end

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

        for at=adapting_temperatures
            PMl=PM;
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
                centered_stimulus = absolute_target - 34 - alpha;
                weight = inv_logit(centered_stimulus.*100);
                theta = (1-weight) * 0.5 + (weight-lambda) .* normcdf(beta .* centered_stimulus);
                
                choice_accuracy=binornd(1,theta);
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                absolute.dataset(n).participant(p).condition(find(adapting_temperatures==at)).baseline_temperature=baseline_temperature;
                absolute.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_adapting_temperature=at;
                
                absolute.dataset(n).participant(p).condition(find(adapting_temperatures==at)).trial(t)=t;
                absolute.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_target_temperature(t)=absolute_target;
                absolute.dataset(n).participant(p).condition(find(adapting_temperatures==at)).choice_accuracy(t)=choice_accuracy;
            end
        end
    end
end
%% Reformat and save the absolute coding simulations as csv
rows={};
row_idx=1;
for n=1:n_datasets
    
    mla=absolute.dataset(n).mu_log_alpha;
    mlb=absolute.dataset(n).mu_log_beta;
    mll=absolute.dataset(n).mu_hlogit_lambda;
    
    tla=absolute.dataset(n).tau_log_alpha;
    tlb=absolute.dataset(n).tau_log_beta;
    tll=absolute.dataset(n).tau_hlogit_lambda; 
        
    for p=1:n_participant
        alpha=absolute.dataset(n).participant(p).alpha;
        beta=absolute.dataset(n).participant(p).beta;
        lambda=absolute.dataset(n).participant(p).lambda;        

        for at_idx=1:length(adapting_temperatures)
            at=adapting_temperatures(at_idx);
            
            for t=1:n_trial
                tt=absolute.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_target_temperature(t);
                ca=absolute.dataset(n).participant(p).condition(find(adapting_temperatures==at)).choice_accuracy(t);
                    
                rows(row_idx,:)={...
                    n,'a',...
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
absolute_table = cell2table(rows, ...
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
writetable(absolute_table, 'absolute_model_data.csv');
%% Simulate trial data for the relative model
for n=1:n_datasets
    fprintf('Simulating relative coding dataset %i out of %i \n',n,n_datasets)
    for p=1:n_participant
        fprintf('Participant %i out of %i \n',p,n_participant)
        alpha=relative.dataset(n).participant(p).alpha;
        beta=relative.dataset(n).participant(p).beta;
        lambda=relative.dataset(n).participant(p).lambda;        

        for at=adapting_temperatures
            PMl=PM;
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
                
                centered_stimulus = relative_target - alpha;
                weight = inv_logit(centered_stimulus.*100);
                theta = (1-weight) * 0.5 + (weight-lambda) .* normcdf(beta .* centered_stimulus);
                
                choice_accuracy=binornd(1,theta);
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                relative.dataset(n).participant(p).condition(find(adapting_temperatures==at)).baseline_temperature=baseline_temperature;
                relative.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_adapting_temperature=at;
                
                relative.dataset(n).participant(p).condition(find(adapting_temperatures==at)).trial(t)=t;
                relative.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_target_temperature(t)=absolute_target;
                relative.dataset(n).participant(p).condition(find(adapting_temperatures==at)).choice_accuracy(t)=choice_accuracy;
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

        for at_idx=1:length(adapting_temperatures)
            at=adapting_temperatures(at_idx);
            
            for t=1:n_trial
                tt=relative.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_target_temperature(t);
                ca=relative.dataset(n).participant(p).condition(find(adapting_temperatures==at)).choice_accuracy(t);
                    
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
writetable(relative_table, 'relative_model_data.csv');
%% Simulate trial data for the mixed model
for n=1:n_datasets
    fprintf('Simulating mixed coding dataset %i out of %i \n',n,n_datasets)
    for p=1:n_participant
        fprintf('Participant %i out of %i \n',p,n_participant)
        alpha=mixed.dataset(n).participant(p).alpha;
        sigma=mixed.dataset(n).participant(p).sigma;
        lambda=mixed.dataset(n).participant(p).lambda;        

        for at=adapting_temperatures
            PMl=PM;
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
                
                centered_stimulus = relative_target - alpha;
                beta = 3/(sigma+2-(at-32));
                weight = inv_logit(centered_stimulus.*100);
                theta = (1-weight) * 0.5 + (weight-lambda) .* normcdf(beta .* centered_stimulus);
                
                choice_accuracy=binornd(1,theta);
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                mixed.dataset(n).participant(p).condition(find(adapting_temperatures==at)).baseline_temperature=baseline_temperature;
                mixed.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_adapting_temperature=at;
                
                mixed.dataset(n).participant(p).condition(find(adapting_temperatures==at)).trial(t)=t;
                mixed.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_target_temperature(t)=absolute_target;
                mixed.dataset(n).participant(p).condition(find(adapting_temperatures==at)).choice_accuracy(t)=choice_accuracy;
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

        for at_idx=1:length(adapting_temperatures)
            at=adapting_temperatures(at_idx);
            
            for t=1:n_trial
                tt=mixed.dataset(n).participant(p).condition(find(adapting_temperatures==at)).absolute_target_temperature(t);
                ca=mixed.dataset(n).participant(p).condition(find(adapting_temperatures==at)).choice_accuracy(t);
                    
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
writetable(mixed_table, 'mixed_model_data.csv');


%%
function y=inv_logit(x) 
    y=1/(1+exp(-x));
end


