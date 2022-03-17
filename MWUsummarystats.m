function [MWUsummary] = MWUsummarystats(MWU_across_all_boots,pvalue_desired,significant_proportion,uniqueindices,Controlsims,Experiment)
%For the matrix MWU_across_all_boots where the rows are genes or ig regions
%and the columns are the number of bootstrap controls that were run on the
%data.  MWU summary is a mX4 matrix where col-1 is the proportion below pval
%column 2 is signicant yes or no where yes=1 and no =0
%column 3 is the average pvalue
%column 4 is the stdev pvalue
[m,n]=size(MWU_across_all_boots);
MWUsummary=zeros(m,6);
for i=1:m
    row=MWU_across_all_boots(i,:);
    instances_below=row(row<=pvalue_desired);
    if isempty(instances_below)~=1;
        proportion=length(instances_below)/n;
        MWUsummary(i,1)=proportion;
        if proportion >= significant_proportion;
            MWUsummary(i,2)=1;
        end
    end
end

%% This section deals with NaN from loci with greater than 10 insertions each and whose differences in reads is indistinguishable (would otherwise by p-val = 1)

for i=1:m;
    norm_MWU=zeros(1,size(MWU_across_all_boots,2));
    for j=1:size(MWU_across_all_boots,2);
        if isnan(MWU_across_all_boots(i,j))==1;
            norm_MWU(1,j)=1;
        else
            norm_MWU(1,j)=MWU_across_all_boots(i,j);
        end
        MWUsummary(i,3)=mean(norm_MWU);
        MWUsummary(i,4)=std(norm_MWU);
    end
end


%% This section gets the average count ratios (and StdDev) between all controlsims and experiment reads.
% The numbers reflect the comparison of Experiment vs Control, so large
% numbers = enriched loci, and low numbers = condtionally essential loci

ctrl_sum=0;
exp_sum=0;
count_ratio=zeros(length(uniqueindices),size(Controlsims,2));

for y=1:size(Controlsims,2); %iterates through the number of columns (simulations) as in controlsims
    for x=1:length(uniqueindices); %goes locus by locus in 1 column
        for p=uniqueindices(x,1):uniqueindices(x,2); %takes the first TAsite and last TA site in the locus
            ctrl_sum=Controlsims(p,y)+ctrl_sum; %replaces the sum as you add in more reads in the locus
            exp_sum=Experiment(p)+exp_sum;
            ratio=exp_sum/ctrl_sum; %calculates the ratio between the samples
        end
        if ratio > 0 %if the ratio for the locus is valid (control and Exp both have reads)
            count_ratio(x,y)=ratio; % deposit the ratio
        end
        if isnan(exp_sum/ctrl_sum)==1; % if the ratio is NaN (0 reads in both control or Exp--e.g., essential locus)
            count_ratio(x,y)==0; % deposit a 0 in the count_ratio matrix
        end
        if isinf(ratio)==1; % if the ratio if Inf (reads in the Exp, but none in the Control)
            count_ratio(x,y)=exp_sum; % deposit the control reads as the ratio (effectively, pretend there is 1 read in exp-->for enriched loci)
        end
        if ratio == 0; %if Exp has no reads, but control does
            count_ratio(x,y)=1/ctrl_sum; % set the ratio as 1/control reads (pretend exp has 1 read).  This is for conditionally essential loci.
        end
    ctrl_sum=0;
    exp_sum=0;
    ratio=0;
    end
end

for z=1:length(count_ratio); % for every locus, calculate the mean fold change and std dev of fold change and report in MWUsummary
    MWUsummary(z,5) = mean(count_ratio(z,:));
    MWUsummary(z,6) = std(count_ratio(z,:));
end
    


    
