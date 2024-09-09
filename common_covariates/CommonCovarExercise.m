% This script considers alternate specifications for the inclusion of
% covariates from W into X

clear all;

% Set sample size

n = 1200;
reps = 100;
seeds = 1:reps;

out_spec1_c0=struct('postsamp',{},'badreps',[]);
out_spec2_c3=struct('postsamp',{},'badreps',[]);

%% Spec 1

% Spec 1 common = 0
parfor (r = 1:reps)

    common = 0;
    %M = strcat('spec1_',num2str(common),'common');
    simdata = GenDataStructureSim2a_CommonCovar_Spec1(seeds(r),n,common);
    try
        postsamp = CollGibbsAlgo_ModelCompRev(simdata.W,simdata.X,simdata.Y,n,simdata.truevals); 
        out_spec1_c0(r).postsamp = postsamp;
    catch ME
        out_spec1_c0(r).badreps = r;
    end
end

% Spec 1 common = 3
parfor (r = 1:reps)

    common = 3;
    simdata = GenDataStructureSim2a_CommonCovar_Spec1(seeds(r),n,common);
    try
        postsamp = CollGibbsAlgo_ModelCompRev(simdata.W,simdata.X,simdata.Y,n,simdata.truevals); 
        out_spec1_c3(r).postsamp = postsamp;
    catch ME
        out_spec1_c3(r).badreps = r;
    end
end

%% Spec 2

% Spec 2 common = 0
parfor (r = 1:reps)

    common = 0;
    %M = strcat('spec1_',num2str(common),'common');
    simdata = GenDataStructureSim2a_CommonCovar_Spec2(seeds(r),n,common);
    try
        postsamp = CollGibbsAlgo_ModelCompRev(simdata.W,simdata.X,simdata.Y,n,simdata.truevals); 
        out_spec2_c0(r).postsamp = postsamp;
    catch ME
        out_spec2_c0(r).badreps = r;
    end
end

% Spec 2 common = 3
parfor (r = 1:reps)

    common = 3;
    simdata = GenDataStructureSim2a_CommonCovar_Spec2(seeds(r),n,common);
    try
        postsamp = CollGibbsAlgo_ModelCompRev(simdata.W,simdata.X,simdata.Y,n,simdata.truevals); 
        out_spec2_c3(r).postsamp = postsamp;
    catch ME
        out_spec2_c3(r).badreps = r;
    end
end

save('common_covars.mat');

%% Report the results

table_spec1_c0 = summarize_commoncovars(out_spec1_c0,4,4,reps);
table_spec1_c3 = summarize_commoncovars(out_spec1_c3,4,4,reps);
table_spec2_c0 = summarize_commoncovars(out_spec2_c0,4,4,reps);
table_spec2_c3 = summarize_commoncovars(out_spec2_c3,4,4,reps);



