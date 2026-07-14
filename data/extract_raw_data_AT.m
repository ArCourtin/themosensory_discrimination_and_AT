% Script used to extract the raw stimulus-response pairs for the AT part of
% the experiment.
% Author: Arthur S. Courtin
clear all

participants_to_extract=[1,3,4,8,12,22,25,26,28,29,30,32,34,36,38,39,40,41,51,56,61];

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

x_r=[];
y_r=[];
t_r=[];
p_r=[];
a_r=[];
c_r=[];
b_r=[];
i_r=[];
bf_r=[];

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
                    load(fullfile(lt(end).folder,lt(end).name));
                    baseline=params.baseline;
                    for tdx=1:length(lt)
                        if contains(lt(tdx).name,'CDT')
                            disp(lt(tdx).name);
                            adapting=str2num(lt(tdx).name(4:end));
                            lc=dir([lt(tdx).folder,'\',lt(tdx).name]);
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
                                    
                                    baseline_tol=0.2;
                                    baseline_flag=zeros(1,l);
                                    for tdx=1:l
                                        temp_trial=PM.temperature{tdx};
                                        temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                        temp_baseline_mean=mean(temp_baseline);
                                        if any(any(abs(temp_baseline-temp_baseline_mean)>baseline_tol))
                                            baseline_flag(tdx)=1;
                                        end
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

                                    baseline_tol=0.2;
                                    baseline_flag=zeros(1,l);
                                    for tdx=1:l
                                        temp_trial=PM.temperature{tdx};
                                        temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                        temp_baseline_mean=mean(temp_baseline);
                                        if any(any(abs(temp_baseline-temp_baseline_mean)>baseline_tol))
                                            baseline_flag(tdx)=1;
                                        end
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
                                end
                            end
                        elseif contains(lt(tdx).name,'WDT')
                            disp(lt(tdx).name);
                            adapting=str2num(lt(tdx).name(4:end));
                            lc=dir([lt(tdx).folder,'\',lt(tdx).name]);
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
                                    
                                    baseline_tol=0.2;
                                    baseline_flag=zeros(1,l);
                                    for tdx=1:l
                                        temp_trial=PM.temperature{tdx};
                                        temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                        temp_baseline_mean=mean(temp_baseline);
                                        if any(any(abs(temp_baseline-temp_baseline_mean)>baseline_tol))
                                            baseline_flag(tdx)=1;
                                        end
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

                                    baseline_tol=0.2;
                                    baseline_flag=zeros(1,l);
                                    for tdx=1:l
                                        temp_trial=PM.temperature{tdx};
                                        temp_baseline=temp_trial(temp_trial(:,1)<0,2:6);
                                        temp_baseline_mean=mean(temp_baseline);
                                        if any(any(abs(temp_baseline-temp_baseline_mean)>baseline_tol))
                                            baseline_flag(tdx)=1;
                                        end
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

fc=table(p_d',b_d',a_d',t_d',i_d',x_d',y_d',ai_d',az_d',rt_d',rep_d',bf_d','VariableNames',{'participant','baseline','adapting','task','trial','temperature','accuracy','active_interval','active_zone','rt','repetition','baseline_flag'});
r=table(p_r',b_r',a_r',t_r',i_r',x_r',y_r',c_r',bf_r','VariableNames',{'participant','baseline','adapting','task','trial','temperature','rating','confirmed','baseline_flag'});

writetable(fc,'d_at_2ifc.csv')
writetable(r,'d_at_ratings.csv')