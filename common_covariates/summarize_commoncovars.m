function evaltab = summarize_commoncovars(out,p,k,reps)


posteval_1sd = zeros(14,reps);
posteval_2sd = zeros(14,reps);
postmeans = zeros(14,reps);
mse = zeros(14,reps);

for r = 1:reps

    % Spec 1 common = 0

    if isempty(out(r).badreps)

        trueVal = out(r).postsamp.Trueval;
        %p = 4; k = 4;
        trueValFlip = zeros(length(trueVal),1);
        trueValFlip(1:p) = trueVal(p+1:2*p); trueValFlip(p+1:2*p) = trueVal(1:p); % flip across the latent classes
        trueValFlip(2*p+1:2*p+k) = -1*trueVal(2*p+1:2*p+k); trueValFlip(2*p+k+1:end) = trueVal(2*p+k+1:end);

        postMean = out(r).postsamp.postMean;
        postSD = out(r).postsamp.postSD;

        cmtrue = trueVal(2*p+1:2*p+k); cmmean = postMean(2*p+1:2*p+k); cmsd = postSD(2*p+1:2*p+k); flipped = zeros(length(cmtrue),1);
        flipped(-cmtrue > -2*cmsd+cmmean & -cmtrue < 2*cmsd+cmmean) =  1;

        flipFlag = 0;
        if sum(flipped) == k
            flipFlag = 1;
        end

        if flipFlag == 0
            posteval_1sd(trueVal > -postSD+postMean & trueVal < postSD+postMean,r) = 1;
            posteval_2sd(trueVal > -2*postSD+postMean & trueVal < 2*postSD+postMean,r) = 1;
            postmeans(:,r) = postMean;
            mse(:,r) = (postmeans(:,r)-trueVal).^2;

        elseif flipFlag == 1
            posteval_1sd(trueValFlip > -postSD+postMean & trueValFlip < postSD+postMean,r) = 1;
            posteval_2sd(trueValFlip > -2*postSD+postMean & trueValFlip < 2*postSD+postMean,r) = 1;
            postmeans(:,r) = [postMean(p+1:2*p);postMean(1:p);-postMean(2*p+1:2*p+k);...
                postMean(2*p+k+2);postMean(2*p+k+1)];
            mse(:,r) = (postmeans(:,r)-trueValFlip).^2;
        end
    else
        posteval_1sd(:,r) = NaN;
        posteval_2sd(:,r) = NaN;
        postmeans(:,r) = NaN;
        mse(:,r) = NaN;

    end

end
varNames = deal(cell(1,6));
varNames(1) = cellstr('parameter');
varNames(2) = cellstr('true value');
varNames(3) = cellstr('posterior mean');
varNames(4) = cellstr('1SD CI cov.');
varNames(5) = cellstr('2SD CI cov.');
varNames(6) = cellstr('mse');
paramvec = {'beta_11','beta_21','beta_31','beta_41',...
    'beta_12','beta_22','beta_32','beta_42',...
    'alpha_1','alpha_2','alpha_3','alpha_4','sig21','sig22'};
evalcell = [paramvec' num2cell(trueVal) num2cell(mean(postmeans,2,'omitnan')) num2cell(mean(posteval_1sd,2,'omitnan'))...
    num2cell(mean(posteval_2sd,2,'omitnan')) num2cell(sum(mse,2,'omitnan'))];
evaltab = cell2table(evalcell);
evaltab.Properties.VariableNames = varNames;
end