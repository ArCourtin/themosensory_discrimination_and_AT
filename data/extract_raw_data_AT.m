% Script used to extract the raw stimulus-response pairs for the AT part of
% the experiment.
% Author: Arthur S. Courtin
clear all

participants_to_extract=[1,3,4,8,12,22,25,26,28,29,30,32,34,36,38,39,40,41,42,51,56,57,61];

x_d=[];
y_d=[];
t_d=[];
p_d=[];
a_d=[];
b_d=[];
i_d=[];
ai_d=[];
az_d=[];
rt_d=[];
rep_d=[];
bf_d=[];
tf_d=[];
df_d=[];
rp_d=[];

x_r=[];
y_r=[];
t_r=[];
p_r=[];
a_r=[];
c_r=[];
b_r=[];
i_r=[];
bf_r=[];
tf_r=[];
df_r=[];
rt_r=[];

lr=dir();
for idx=1:length(lr)
    if char(lr(idx).name(1))=='s'
        participant=str2num(lr(idx).name(5:8));
        disp(lr(idx).name);
        if sum(participant==participants_to_extract)
            lp=dir([lr(idx).folder,'\',lr(idx).name]);
            ld=dir([lp(3).folder,'\',lp(3).name]);
            for ddx=1:length(ld)
                if strcmp(ld(ddx).name,{'D_AT'})
                    lt=dir([ld(ddx).folder,'\',ld(ddx).name]);
                    params=[];
                    baseline=[];
                    adaptation_choice=[];
                    load(fullfile(lt(end).folder,lt(end).name));
                    baseline=params.baseline;
                    adaptation_choice=adaptation_choice.adaptation_choice;
                    for tdx=1:length(lt)
                        if contains(lt(tdx).name,'CDT')
                            disp(lt(tdx).name);
                            adapting=str2num(lt(tdx).name(4:end));
                            lc=dir([lt(tdx).folder,'\',lt(tdx).name]);
                            if(adaptation_choice(find(adapting==params.new_baseline))==1)
                                for cdx=1:length(lc)
                                    if contains(lc(cdx).name,'2ifc')&&contains(lc(cdx).name,'.mat')
                                        PM=[];
                                        load([lc(cdx).folder,'\',lc(cdx).name]);
                                        y=PM.response;
                                        l=length(y);
                                        x=PM.x(1:l);
                                        ai=PM.active_interval(1:l);
                                        az=PM.active_zones(1:l);
                                        rt=PM.rt(1:l);
                                        rep=PM.rep(1:l);
                                        
                                        tolerance=0.2;
                                        baseline_flag=zeros(1,l);
                                        target_flag=zeros(1,l);
                                        deviation_flag=nan(1,l);
                                        recorded_temp=nan(1,l);
                                        for tdx=1:l
                                            temp_trial=PM.temperature{tdx};
                                            
                                            temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                            baseline_flag(tdx)=any(any(abs(temp_baseline-adapting)>tolerance));
                                          
                                            target=adapting-PM.x(tdx);
                                            temp_plateau=temp_trial((temp_trial(:,1)<2)&(temp_trial(:,1)>1.5),2:6);
                                            temp_plateau_mean=mean(temp_plateau);
                                            if PM.active_zones(tdx)==1
                                                zdx=[4 5];
                                            else
                                                zdx=[1 2];
                                            end
                                            target_flag(tdx)=any(abs(temp_plateau_mean(zdx)-target)>tolerance);
                                            deviation_flag(tdx)=any(any(abs(temp_plateau_mean-temp_plateau)>tolerance));
                                            recorded_temp(tdx)=mean(temp_plateau_mean(zdx));
                                        end
    
                                        t=zeros(1,l);
                                        p=repmat(participant,1,l);
                                        a=repmat(adapting,1,l);
                                        b=repmat(baseline,1,l);
                                        i=1:l;
                                        
                                        x_d=[x_d x];
                                        y_d=[y_d y];
                                        ai_d=[ai_d ai];
                                        az_d=[az_d az];
                                        rt_d=[rt_d rt];
                                        rep_d=[rep_d rep];
                                        bf_d=[bf_d baseline_flag];
                                        tf_d=[tf_d target_flag];
                                        df_d=[df_d deviation_flag];
                                        rp_d=[rp_d recorded_temp];
                                        t_d=[t_d t];
                                        p_d=[p_d p];                        
                                        a_d=[a_d a];
                                        b_d=[b_d b];
                                        i_d=[i_d i];
                                       
                                    elseif contains(lc(cdx).name,'ratings')&&contains(lc(cdx).name,'.mat')
                                        ratings=[];
                
                                        load([lc(cdx).folder,'\',lc(cdx).name]);
                
                                        y=ratings.rating;
                                        l=length(y);
                                        x=ratings.target;
                                        c=ratings.confirmed;
                                        t=zeros(1,l);
                                        p=repmat(participant,1,l);
                                        a=repmat(adapting,1,l);
                                        b=repmat(baseline,1,l);
                                        i=1:l;
    
                                        tolerance=0.2;
                                        baseline_flag=zeros(1,l);
                                        target_flag=zeros(1,l);
                                        deviation_flag=nan(1,l);
                                        recorded_temp=nan(1,l);
                                        for tdx=1:l
                                            temp_trial=ratings.temperature{tdx};
                                            
                                            temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                            baseline_flag(tdx)=any(any(abs(temp_baseline-adapting)>tolerance));
                                          
                                            target=ratings.target(tdx);
                                            temp_plateau=temp_trial((temp_trial(:,1)<2)&(temp_trial(:,1)>1.5),2:6);
                                            temp_plateau_mean=mean(temp_plateau);
                                            if mod(tdx,2)
                                                zdx=[4 5];
                                            else
                                                zdx=[1 2];
                                            end
                                            target_flag(tdx)=any(abs(temp_plateau_mean(zdx)-target)>tolerance);
                                            deviation_flag(tdx)=any(any(abs(temp_plateau_mean-temp_plateau)>tolerance));
                                            recorded_temp(tdx)=mean(temp_plateau_mean(zdx));
                                        end              

                                        x_r=[x_r x];
                                        y_r=[y_r y];
                                        t_r=[t_r t];
                                        p_r=[p_r p];                        
                                        a_r=[a_r a];
                                        c_r=[c_r c];
                                        i_r=[i_r i];
                                        b_r=[b_r b];
                                        bf_r=[bf_r baseline_flag];
                                        tf_r=[tf_r target_flag];
                                        df_r=[df_r deviation_flag];
                                        rt_r=[rt_r recorded_temp];
                                    end
                                end
                            end
                        elseif contains(lt(tdx).name,'WDT')
                            disp(lt(tdx).name);
                            adapting=str2num(lt(tdx).name(4:end));
                            lc=dir([lt(tdx).folder,'\',lt(tdx).name]);
                            if(adaptation_choice(find(adapting==params.new_baseline))==1)
                                for cdx=1:length(lc)
                                    if contains(lc(cdx).name,'2ifc')&&contains(lc(cdx).name,'.mat')
                                        PM=[];
                                        load([lc(cdx).folder,'\',lc(cdx).name]);
                                        y=PM.response;
                                        l=length(y);
                                        x=PM.x(1:l);
                                        ai=PM.active_interval(1:l);
                                        az=PM.active_zones(1:l);
                                        rt=PM.rt(1:l);
                                        rep=PM.rep(1:l);
                                        
                                        tolerance=0.2;
                                        baseline_flag=zeros(1,l);
                                        target_flag=zeros(1,l);
                                        deviation_flag=nan(1,l);
                                        recorded_temp=nan(1,l);
                                        for tdx=1:l
                                            temp_trial=PM.temperature{tdx};
                                            
                                            temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                            baseline_flag(tdx)=any(any(abs(temp_baseline-adapting)>tolerance));
                                          
                                            target=adapting+PM.x(tdx);
                                            temp_plateau=temp_trial((temp_trial(:,1)<2)&(temp_trial(:,1)>1.5),2:6);
                                            temp_plateau_mean=mean(temp_plateau);
                                            if PM.active_zones(tdx)==1
                                                zdx=[4 5];
                                            else
                                                zdx=[1 2];
                                            end
                                            target_flag(tdx)=any(abs(temp_plateau_mean(zdx)-target)>tolerance);
                                            deviation_flag(tdx)=any(any(abs(temp_plateau_mean-temp_plateau)>tolerance));
                                            recorded_temp(tdx)=mean(temp_plateau_mean(zdx));
                                        end

                                        t=ones(1,l);
                                        p=repmat(participant,1,l);
                                        a=repmat(adapting,1,l);
                                        b=repmat(baseline,1,l);
                                        i=1:l;
                                        
                                        x_d=[x_d x];
                                        y_d=[y_d y];
                                        t_d=[t_d t];
                                        ai_d=[ai_d ai];
                                        az_d=[az_d az];
                                        rt_d=[rt_d rt];
                                        rep_d=[rep_d rep];
                                        bf_d=[bf_d baseline_flag];
                                        tf_d=[tf_d target_flag];
                                        df_d=[df_d deviation_flag];
                                        rp_d=[rp_d recorded_temp];
                                        p_d=[p_d p];                        
                                        a_d=[a_d a];
                                        b_d=[b_d b];
                                        i_d=[i_d i];
                                       
                                    elseif contains(lc(cdx).name,'ratings')&&contains(lc(cdx).name,'.mat')
                                        ratings=[];
                
                                        load([lc(cdx).folder,'\',lc(cdx).name]);
                
                                        y=ratings.rating;
                                        l=length(y);
                                        x=ratings.target;
                                        c=ratings.confirmed;
                                        t=ones(1,l);
                                        p=repmat(participant,1,l);
                                        a=repmat(adapting,1,l);
                                        b=repmat(baseline,1,l);
                                        i=1:l;
    
                                        tolerance=0.2;
                                        baseline_flag=zeros(1,l);
                                        target_flag=zeros(1,l);
                                        deviation_flag=nan(1,l);
                                        recorded_temp=nan(1,l);
                                        for tdx=1:l
                                            temp_trial=ratings.temperature{tdx};
                                            
                                            temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                            baseline_flag(tdx)=any(any(abs(temp_baseline-adapting)>tolerance));
                                          
                                            target=ratings.target(tdx);
                                            temp_plateau=temp_trial((temp_trial(:,1)<2)&(temp_trial(:,1)>1.5),2:6);
                                            temp_plateau_mean=mean(temp_plateau);
                                            if mod(tdx,2)
                                                zdx=[4 5];
                                            else
                                                zdx=[1 2];
                                            end
                                            target_flag(tdx)=any(abs(temp_plateau_mean(zdx)-target)>tolerance);
                                            deviation_flag(tdx)=any(any(abs(temp_plateau_mean-temp_plateau)>tolerance));
                                            recorded_temp(tdx)=mean(temp_plateau_mean(zdx));
                                        end
    
                                        x_r=[x_r x];
                                        y_r=[y_r y];
                                        t_r=[t_r t];
                                        p_r=[p_r p];                        
                                        a_r=[a_r a];
                                        c_r=[c_r c];
                                        i_r=[i_r i];
                                        b_r=[b_r b];
                                        bf_r=[bf_r baseline_flag];
                                        tf_r=[tf_r target_flag];
                                        df_r=[df_r deviation_flag];
                                        rt_r=[rt_r recorded_temp];
                                    end
                                end
                            end                        
                        end
                    end
                else
                end
            end
        end
    end
end

fc=table(p_d',b_d',a_d',t_d',i_d',x_d',y_d',ai_d',az_d',rt_d',rep_d',bf_d',tf_d',df_d',rp_d','VariableNames',{'participant','baseline','adapting','task','trial','temperature','accuracy','active_interval','active_zone','rt','repetition','baseline_flag','target_flag','deviation_flag','recorded_temperature'});
r=table(p_r',b_r',a_r',t_r',i_r',x_r',y_r',c_r',bf_r',tf_r',df_r',rt_r','VariableNames',{'participant','baseline','adapting','task','trial','temperature','rating','confirmed','baseline_flag','target_flag','deviation_flag','recorded_temperature'});

writetable(fc,'d_at_2ifc_af.csv')
writetable(r,'d_at_ratings_af.csv')