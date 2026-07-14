%Script to simulate data from the "mechanistic" discrimination models
%Licence: MIT
%Author: Arthur S. Courtin
%Edited with the assistance of Claude Code (Anthropic).

%% General settings
rng(12345)
addpath(genpath(fullfile(fileparts(cd), 'Palamedes')));

n_datasets = 30;
n_participant = 20;
n_trial = 30;
n_rating_trial = 10;

baseline_temperature = 32;
adapting_temperatures = baseline_temperature + (-2:2);

plot_pp=0;

%% Generate group and participant parameters for the absolute model
close all
for n=1:n_datasets
    rho(1:n_participant)=-10;
    while any(rho<-2)
        m_r=-10;
        while m_r<-2
            m_r=normrnd(0,2);
        end
        absolute.dataset(n).mu_rho=m_r;
        absolute.dataset(n).mu_log_alpha=normrnd(-2,1);
        absolute.dataset(n).mu_log_beta=normrnd(0,1);
        absolute.dataset(n).mu_hlogit_lambda=normrnd(-4,1);
        absolute.dataset(n).mu_kappa=normrnd(0,1);

        absolute.dataset(n).mu_intercept=normrnd(-2,1);
        absolute.dataset(n).mu_log_slope=normrnd(-2,1);
        absolute.dataset(n).mu_log_lb=normrnd(2,.5);
        absolute.dataset(n).mu_log_ub=normrnd(2,.5);
        absolute.dataset(n).mu_log_eta=normrnd(3,1);

        absolute.dataset(n).tau_rho=abs(normrnd(0,2));
        absolute.dataset(n).tau_log_alpha=abs(normrnd(0,1));
        absolute.dataset(n).tau_log_beta=abs(normrnd(0,1));
        absolute.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));
        absolute.dataset(n).tau_kappa=abs(normrnd(0,1));

        absolute.dataset(n).tau_intercept=abs(normrnd(-2,1));
        absolute.dataset(n).tau_log_slope=abs(normrnd(-2,1));
        absolute.dataset(n).tau_log_lb=abs(normrnd(2,.5));
        absolute.dataset(n).tau_log_ub=abs(normrnd(2,.5));
        absolute.dataset(n).tau_log_eta=abs(normrnd(3,1));

        for p=1:n_participant
            rho(p)=normrnd(absolute.dataset(n).mu_rho,absolute.dataset(n).tau_rho);
            absolute.dataset(n).participant(p).rho=rho(p);
            absolute.dataset(n).participant(p).alpha=exp(normrnd(absolute.dataset(n).mu_log_alpha,absolute.dataset(n).tau_log_alpha));
            absolute.dataset(n).participant(p).beta=exp(normrnd(absolute.dataset(n).mu_log_beta,absolute.dataset(n).tau_log_beta));
            absolute.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(absolute.dataset(n).mu_hlogit_lambda,absolute.dataset(n).tau_hlogit_lambda)));
            absolute.dataset(n).participant(p).kappa=normrnd(absolute.dataset(n).mu_kappa,absolute.dataset(n).tau_kappa);

            absolute.dataset(n).participant(p).intercept=normrnd(absolute.dataset(n).mu_intercept,absolute.dataset(n).tau_intercept);
            absolute.dataset(n).participant(p).slope=exp(normrnd(absolute.dataset(n).mu_log_slope,absolute.dataset(n).tau_log_slope));
            absolute.dataset(n).participant(p).lb=exp(normrnd(absolute.dataset(n).mu_log_lb,absolute.dataset(n).tau_log_lb));
            absolute.dataset(n).participant(p).ub=exp(normrnd(absolute.dataset(n).mu_log_ub,absolute.dataset(n).tau_log_ub));
            absolute.dataset(n).participant(p).eta=exp(normrnd(absolute.dataset(n).mu_log_eta,absolute.dataset(n).tau_log_eta));
        end
    end
    if plot_pp
        x=30:.1:40;
        d=-10:.1:10;
        figure
        hold on
        for p=1:n_participant
            alpha=absolute.dataset(n).participant(p).alpha;
            rho=absolute.dataset(n).participant(p).rho;
            beta=absolute.dataset(n).participant(p).beta;
            lambda=absolute.dataset(n).participant(p).lambda;
            kappa=absolute.dataset(n).participant(p).kappa;
    
            for at=adapting_temperatures
                % Response-coded 2IFC: d is the signed target deviation, interval_sign flags which
                % interval held the deviating stimulus and kappa is the interval bias (matches the
                % simulated model below and plot_priors.R). y = P(chose 2nd interval).
                interval_sign=sign(d);
                target=at+abs(d);
                x_c=target-baseline_temperature-rho;
                stim_rep=x_c./(1+exp(-100*x_c));
                adapt_stim_rep=stim_rep./(1+exp(-100*(stim_rep-(at-baseline_temperature+alpha-rho))));
                theta=lambda+(1-2*lambda)*normcdf(interval_sign.*beta.*adapt_stim_rep-kappa);
                plot(d,theta)
            end
        end
        hold off
        
        figure
        hold on
        for p=1:n_participant
            intercept=absolute.dataset(n).participant(p).intercept;
            slope=absolute.dataset(n).participant(p).slope;
          
            for at=adapting_temperatures
                x_c=x-baseline_temperature;
                theta=1./(1+exp(-(intercept+slope*x_c)));
                plot(x,theta)
            end
        end
        hold off
    end
end
%% Generate group and participant parameters for the relative model
close all
for n=1:n_datasets
    relative.dataset(n).mu_log_alpha=normrnd(-2,1);
    relative.dataset(n).mu_log_beta=normrnd(0,1);
    relative.dataset(n).mu_hlogit_lambda=normrnd(-4,1);
    relative.dataset(n).mu_kappa=normrnd(0,1);
    
    relative.dataset(n).mu_intercept=normrnd(-2,1);
    relative.dataset(n).mu_log_slope=normrnd(-2,1);
    relative.dataset(n).mu_log_lb=normrnd(2,.5);
    relative.dataset(n).mu_log_ub=normrnd(2,.5);
    relative.dataset(n).mu_log_eta=normrnd(3,1);
    
    relative.dataset(n).tau_log_alpha=abs(normrnd(0,1));
    relative.dataset(n).tau_log_beta=abs(normrnd(0,1));
    relative.dataset(n).tau_hlogit_lambda=abs(normrnd(0,1));
    relative.dataset(n).tau_kappa=abs(normrnd(0,1));
    
    relative.dataset(n).tau_intercept=abs(normrnd(-2,1));
    relative.dataset(n).tau_log_slope=abs(normrnd(-2,1));
    relative.dataset(n).tau_log_lb=abs(normrnd(2,.5));
    relative.dataset(n).tau_log_ub=abs(normrnd(2,.5));
    relative.dataset(n).tau_log_eta=abs(normrnd(3,1));
    
    for p=1:n_participant
        relative.dataset(n).participant(p).alpha=exp(normrnd(relative.dataset(n).mu_log_alpha,relative.dataset(n).tau_log_alpha));
        relative.dataset(n).participant(p).beta=exp(normrnd(relative.dataset(n).mu_log_beta,relative.dataset(n).tau_log_beta));
        relative.dataset(n).participant(p).lambda=.5 /(1+exp(-normrnd(relative.dataset(n).mu_hlogit_lambda,relative.dataset(n).tau_hlogit_lambda)));
        relative.dataset(n).participant(p).kappa=normrnd(relative.dataset(n).mu_kappa,relative.dataset(n).tau_kappa);
        
        relative.dataset(n).participant(p).intercept=normrnd(relative.dataset(n).mu_intercept,relative.dataset(n).tau_intercept);
        relative.dataset(n).participant(p).slope=exp(normrnd(relative.dataset(n).mu_log_slope,relative.dataset(n).tau_log_slope));
        relative.dataset(n).participant(p).lb=exp(normrnd(relative.dataset(n).mu_log_lb,relative.dataset(n).tau_log_lb));
        relative.dataset(n).participant(p).ub=exp(normrnd(relative.dataset(n).mu_log_ub,relative.dataset(n).tau_log_ub));
        relative.dataset(n).participant(p).eta=exp(normrnd(relative.dataset(n).mu_log_eta,relative.dataset(n).tau_log_eta));
    end
    
    if plot_pp
        d=-10:.1:10;
        figure
        hold on
        for p=1:n_participant
        alpha=relative.dataset(n).participant(p).alpha;
        beta=relative.dataset(n).participant(p).beta;
        lambda=relative.dataset(n).participant(p).lambda;
        kappa=relative.dataset(n).participant(p).kappa;
    
            for at=adapting_temperatures
                % Response-coded 2IFC: d is the signed target deviation, interval_sign flags which
                % interval held the deviating stimulus and kappa is the interval bias (matches the
                % simulated model below and plot_priors.R). Relative coding is at-invariant, so the
                % five adapting-temperature curves overlap. y = P(chose 2nd interval).
                interval_sign=sign(d);
                x_c=abs(d)-alpha;
                stim_rep=x_c./(1+exp(-100*x_c));
                theta=lambda+(1-2*lambda)*normcdf(interval_sign.*beta.*stim_rep-kappa);
                plot(d,theta)
            end
        end
        hold off
        x=30:.1:40;
        figure
        hold on
        for p=1:n_participant
        intercept=relative.dataset(n).participant(p).intercept;
        slope=relative.dataset(n).participant(p).slope;
            for at=adapting_temperatures
                x_c=x-at;
                theta=1./(1+exp(-(intercept+slope*x_c)));
                plot(x,theta)
            end
        end
        hold off
    end
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
        rho=absolute.dataset(n).participant(p).rho;
        beta=absolute.dataset(n).participant(p).beta;
        lambda=absolute.dataset(n).participant(p).lambda;        
        
        intercept=absolute.dataset(n).participant(p).intercept;
        slope=absolute.dataset(n).participant(p).slope;
        lb=absolute.dataset(n).participant(p).lb;
        ub=absolute.dataset(n).participant(p).ub;
        eta=absolute.dataset(n).participant(p).eta;

        for at=adapting_temperatures
            adx=at-baseline_temperature+3;
            PMl=PM;
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
                stimulus_representation = (absolute_target-baseline_temperature-rho)/(1+exp(-100*(absolute_target-baseline_temperature-rho)));
                adapted_stimulus_representation = stimulus_representation/(1+exp(-100*(stimulus_representation-(at-baseline_temperature+alpha-rho))));
                % Response-coded 2IFC. Which interval holds the deviating stimulus is exactly
                % balanced within each condition: generated up front in blocks of 6 (3 first-interval,
                % 3 second-interval), matching the experiment. The second-interval choice is drawn
                % under the interval bias, and accuracy is derived from it to feed the (unchanged) Psi staircase.
                if t==1
                    active_interval_sequence=[];
                    for b=1:(n_trial/6)
                        block=[1 1 1 2 2 2];
                        active_interval_sequence=[active_interval_sequence block(randperm(6))];
                    end
                end
                active_interval = active_interval_sequence(t);
                interval_sign = 2*(active_interval==2)-1;
                p_second = lambda+(1-2*lambda)*normcdf(interval_sign*beta*adapted_stimulus_representation-absolute.dataset(n).participant(p).kappa);
                chose_second = binornd(1,p_second);
                if active_interval==2
                    choice_accuracy = chose_second;
                else
                    choice_accuracy = 1-chose_second;
                end
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                absolute.dataset(n).participant(p).condition(adx).baseline_temperature=baseline_temperature;
                absolute.dataset(n).participant(p).condition(adx).absolute_adapting_temperature=at;
                
                absolute.dataset(n).participant(p).condition(adx).trial(t)=t;
                absolute.dataset(n).participant(p).condition(adx).absolute_target_temperature(t)=absolute_target;
                absolute.dataset(n).participant(p).condition(adx).choice_accuracy(t)=choice_accuracy;
                absolute.dataset(n).participant(p).condition(adx).active_interval(t)=active_interval;
            end
            
            [target_temperatures] = determine_target_detection_at(PMl,at);
            for t=1:n_rating_trial
                if mod(t,2)
                    tt=target_temperatures(1);
                else
                    tt=target_temperatures(2);
                end
               
                lr=intercept+slope*(tt-baseline_temperature);

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
    
    mr=absolute.dataset(n).mu_rho;
    mlb=absolute.dataset(n).mu_log_beta;
    mla=absolute.dataset(n).mu_log_alpha;
    mll=absolute.dataset(n).mu_hlogit_lambda;
    mk=absolute.dataset(n).mu_kappa;

    tr=absolute.dataset(n).tau_rho;
    tlb=absolute.dataset(n).tau_log_beta;
    tla=absolute.dataset(n).tau_log_alpha;
    tll=absolute.dataset(n).tau_hlogit_lambda;
    tk=absolute.dataset(n).tau_kappa;
        
    for p=1:n_participant
        alpha=absolute.dataset(n).participant(p).alpha;
        beta=absolute.dataset(n).participant(p).beta;
        rho=absolute.dataset(n).participant(p).rho;
        lambda=absolute.dataset(n).participant(p).lambda;
        kappa=absolute.dataset(n).participant(p).kappa;

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_trial
                tt=absolute.dataset(n).participant(p).condition(adx).absolute_target_temperature(t);
                ca=absolute.dataset(n).participant(p).condition(adx).choice_accuracy(t);
                ai=absolute.dataset(n).participant(p).condition(adx).active_interval(t);

                rows(row_idx,:)={...
                    n,'a',...
                    mr,mla,mlb,mll,mk,...
                    tr,tla,tlb,tll,tk,...
                    p,...
                    alpha,rho,beta,lambda,kappa,...
                    baseline_temperature,at,...
                    t,tt,ca,ai...
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
        'mu_rho',...
        'mu_log_alpha',...
        'mu_log_beta',...
        'mu_hlogit_lambda',...
        'mu_kappa',...
        'tau_rho',...
        'tau_log_alpha',...
        'tau_log_beta',...
        'tau_hlogit_lambda',...
        'tau_kappa',...
        'participant', ...
        'alpha', ...
        'rho', ...
        'beta', ...
        'lambda', ...
        'kappa', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'choice_accuracy', ...
        'active_interval' ...
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
            adx=at-baseline_temperature+3;
            PMl=PM;
            
            for t=1:n_trial
                relative_target = PMl.xCurrent;
                absolute_target = relative_target + at;
                stimulus_representation = (relative_target-alpha)/(1+exp(-100*(relative_target-alpha)));
                % Response-coded 2IFC. Which interval holds the deviating stimulus is exactly
                % balanced within each condition: generated up front in blocks of 6 (3 first-interval,
                % 3 second-interval), matching the experiment. The second-interval choice is drawn
                % under the interval bias, and accuracy is derived from it to feed the (unchanged) Psi staircase.
                if t==1
                    active_interval_sequence=[];
                    for b=1:(n_trial/6)
                        block=[1 1 1 2 2 2];
                        active_interval_sequence=[active_interval_sequence block(randperm(6))];
                    end
                end
                active_interval = active_interval_sequence(t);
                interval_sign = 2*(active_interval==2)-1;
                p_second = lambda+(1-2*lambda)*normcdf(interval_sign*beta*stimulus_representation-relative.dataset(n).participant(p).kappa);
                chose_second = binornd(1,p_second);
                if active_interval==2
                    choice_accuracy = chose_second;
                else
                    choice_accuracy = 1-chose_second;
                end
                PMl=PAL_AMPM_updatePM(PMl,choice_accuracy);
                
                relative.dataset(n).participant(p).condition(adx).baseline_temperature=baseline_temperature;
                relative.dataset(n).participant(p).condition(adx).absolute_adapting_temperature=at;
                
                relative.dataset(n).participant(p).condition(adx).trial(t)=t;
                relative.dataset(n).participant(p).condition(adx).absolute_target_temperature(t)=absolute_target;
                relative.dataset(n).participant(p).condition(adx).choice_accuracy(t)=choice_accuracy;
                relative.dataset(n).participant(p).condition(adx).active_interval(t)=active_interval;
            end
            
            [target_temperatures] = determine_target_detection_at(PMl,at);
            for t=1:n_rating_trial
                if mod(t,2)
                    tt=target_temperatures(1);
                else
                    tt=target_temperatures(2);
                end

                lr=intercept+slope*(tt-at);

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
    mk=relative.dataset(n).mu_kappa;

    tla=relative.dataset(n).tau_log_alpha;
    tlb=relative.dataset(n).tau_log_beta;
    tll=relative.dataset(n).tau_hlogit_lambda;
    tk=relative.dataset(n).tau_kappa;
        
    for p=1:n_participant
        alpha=relative.dataset(n).participant(p).alpha;
        beta=relative.dataset(n).participant(p).beta;
        lambda=relative.dataset(n).participant(p).lambda;
        kappa=relative.dataset(n).participant(p).kappa;

        for adx=1:length(adapting_temperatures)
            at=adapting_temperatures(adx);
            
            for t=1:n_trial
                tt=relative.dataset(n).participant(p).condition(adx).absolute_target_temperature(t);
                ca=relative.dataset(n).participant(p).condition(adx).choice_accuracy(t);
                ai=relative.dataset(n).participant(p).condition(adx).active_interval(t);

                rows(row_idx,:)={...
                    n,'r',...
                    mla,mlb,mll,mk,...
                    tla,tlb,tll,tk,...
                    p,...
                    alpha,beta,lambda,kappa,...
                    baseline_temperature,at,...
                    t,tt,ca,ai...
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
        'mu_kappa',...
        'tau_log_alpha',...
        'tau_log_beta',...
        'tau_hlogit_lambda',...
        'tau_kappa',...
        'participant', ...
        'alpha', ...
        'beta', ...
        'lambda', ...
        'kappa', ...
        'recorded_baseline_temperature', ...
        'absolute_adapting_temperature', ...
        'trial', ...
        'absolute_target_temperature', ...
        'choice_accuracy', ...
        'active_interval' ...
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

