function annotate_Akk_TAsites_Artist(TAsites,gbk_file)
%last edit: Feb-2020

%Applies annotations to TA sites (column vector of genome coordinates)
%INPUT
    % TA sites in genome
    % Annotation file (gbff or gbk)
%OUTPUT
    % Add geneID to table and save as new file
    % GeneID can be: CDS (Amuc_####), intergenic (IG_Amuc_####)
        %convention is to name intergenics by the gene they preceed
    
 %geneID / class / orientation / size / TN.position / product / notes
        %geneID = Amuc# (locus_tag); Amuc_M####
        %class = CDS, tRNA, rRNA, ncRNA, pseudo, intergenic, misc
        %orientation = -, +, na
        %size = distance in bp btwn indices
        %TN.position = relative position of transposon
            %relative to orientation, when orientaton is given
            %relative to ascending coordinate when orientation 'na'
        %product = description of predicted product, na
        %[notes = class specific information, e.g. sig-peptide for CDS]
    
%PROCESS
    %For each CDS annotation query table, collect matching indices & data

%IMPROVE
    %Add "misc_binding" feature information = riboswitch (2x in Akk genome)
    
    
%%
%Format inputs  
   %NEEDs work to accommodate different types of tables: eg. -/+ header
%gbk_file='/Users/pmalkus/Documents/Valdivia_Lab/Akkermansia/INseq_stuff/Annotation_stuff/Amuc_CP001071.1.gb.txt';
gbkStruct = genbankread(gbk_file);
f = featureparse(gbkStruct);
%read 1,2 from results table and store as cell array

%load TA sites (assume format from output of "TAfinder.py")
coor = dlmread(TAsites,'',0,1);

%Holders for: geneID, class, orientation, size, rel.position,...
  %...product decription, signal peptide (-/+), notes
len = length(coor);
%annot = cell(length(coor),7);
gene=cell(len,1); class=cell(len,1); orient=cell(len,1);
siz=[len,1];pos=[len,1];
prod=cell(len,1);sigP=cell(len,1);note=cell(len,1);

%Collect all 'gene' info (assign pseudo); define intergenic regions
  %Use to cycle through other 'features' for matching 'product'(or 'note')
sq=f.gene;
indices=[sq.Indices];
%FIX indices from gbk file to match corrected genome (remove G at 1704819)
modi=find(indices < 1704819,1,'last');
indices(modi+1:end)=indices(modi+1:end)-1;
indices=[indices(1:2:end);indices(2:2:end)];%gene coordinates, [start stop]
indices=indices';

loc=sort(indices,2,'ascend');%gene coordinates (location), [low high]  
for i=1:length(coor)%for all coordinates in results table
    c=coor(i);%retrieve coordinate
    gidx=find(loc(:,1)<c,1,'last');%retrieve locus before insertion
    %Intergenic handler (assign all intergenic as orientation = '+')
    if isempty(gidx)%special case: insertion before first CDS
        name1 = {'IG'}; name2 = {sq(1).locus_tag};
        gene(i) = {strjoin([name1 name2],'_')};
        siz(i)=loc(1);%assumes first coordinate =1
        pos(i) = c/siz(i);
        class(i) = {'intergenic'}; orient(i) = {'+'};
    elseif c > loc(gidx,2) %is intergenic?
        if gidx==length(sq) %special case: insertion after last gene
            gene(i) = {'IG_end'};
            siz(i)=f.source.Indices(2)-loc(gidx,2);
            pos(i) = (c-loc(gidx,2))/siz(i);
            class(i) = {'intergenic'}; orient(i) = {'+'};
            continue
        end
        name1 = {'IG'}; name2 = {sq(gidx+1).locus_tag};
        gene(i) = {strjoin([name1 name2],'_')};
        siz(i) = loc(gidx+1,1)-loc(gidx,2);
        pos(i) = (c-loc(gidx,2))/siz(i);
        class(i) = {'intergenic'}; orient(i) = {'+'}; 
    else 
        %for genes: extract & store info
        gene(i) = {sq(gidx).locus_tag};
        siz(i) = loc(gidx,2)-loc(gidx,1);
        pos(i) = (abs(c-indices(gidx,1)))/siz(i);
        if sq(gidx).pseudo==1
            class(i) = {'pseudo'};
        else
            class(i) = {'gene'};
        end
        if indices(gidx,2) > indices(gidx,1)
            orient(i) = {'+'};
        else
            orient(i) = {'-'};
        end
    end
end
siz=siz'; pos=pos';

%Make a table with relevant fields for sub-structure of Features (f)
  %class, locus_tag, product, note --> from f.tRNA, f.CDS, f.rRNA, f.ncRNA 
classID(1:length(f.tRNA),1) = {'tRNA'};
classID(end:end+length(f.CDS),1) = {'CDS'};
classID(end:end+length(f.rRNA),1) = {'rRNA'};
classID(end:end+length(f.ncRNA),1) = {'ncRNA'};
locusID = [{f.tRNA.locus_tag}';{f.CDS.locus_tag}';{f.rRNA.locus_tag}';{f.ncRNA.locus_tag}'];
prodID = [{f.tRNA.product}';{f.CDS.product}';{f.rRNA.product}';{f.ncRNA.note}'];

%Find annotation by locus_tag; add to holder cell arrays
    %product decription, signal peptide (-/+), notes
for j=1:len
    g=gene(j);
    lidx=find(strcmp(locusID,g));
    if isempty(lidx)
        prod(j)={'n/a'};
        continue
    end
    prod(j)=prodID(lidx);
    if strcmp(class(j),{'gene'})==1
        class(j)=classID(lidx);
    end
end

%Signal peptide information
sig=f.sig_peptide;
siglocusID = {sig.locus_tag}';
for j=1:len
    g=gene(j);
    lidx=find(strcmp(siglocusID,g));
    if isempty(lidx)
        sigP(j) = {'-'};
    else
        sigP(j) = {'+'};
        %note(j) = {sig(lidx).note};
    end
end
    
%Add annotations to results table
table = num2cell(coor);
csiz = num2cell(siz);
cpos = num2cell(pos);
out = [table gene class orient csiz cpos prod sigP];%removed 'note'
%Save output as tab-delimited file
passfile = strcat(TAsites(1:end-4),'ID.txt');
%writecell(out,passfile,'Delimiter','tab');%for 2019a

% Convert cell to a table and use first row as variable names
%T = cell2table(out(2:end,:),'VariableNames',out(1,:))
% Write the table to a CSV file
%writetable(T,passfile)

%write cell array to tab-delimted text
fileID = fopen(passfile,'w');
formatSpec = '%d\t %s\t %s\t %s\t %d\t %d\t %s\t %s\n';
[nrows,ncols] = size(out);
for row = 1:nrows
    fprintf(fileID,formatSpec,out{row,:});
end
fclose(fileID);
    
end
