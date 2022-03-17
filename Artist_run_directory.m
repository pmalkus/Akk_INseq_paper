function Artist_run_directory(genome,inputName,filterName)
%last edit: May-13-2020, works

% Acts on all Test Samples in current directory

% Consolidates ARTIST workflow
% Generate normalized Input (mean of ControlSims)

% Parse all inputs and generate .mat file containing all workspace
    % variables needed for conditional gene analysis with ARTIST 
    
% INPUTs
    % 'Genome', prefix for Genome_TAsites and Genome_TAsitesID .txt files
    % 'inputName', filname: sample name for read table from Goodman pipeline
        % coordinates in column-2, total reads in column-5
    % 'filterName', sample name for read table from Goodman pipeline
        % used to filter sample data prior to filling missing TA sites

% OUTPUT
    % .mat files named by inputName_sampleName_filterName
    % contains workspace variables for analysis and MWU stats
        % also TAhits: Tn insertions per locus
        % and NormInput: mean of simulated/normalized counts for each site

        
% Loop through all Sample files in directory, excluding Input and Filter files
    
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove Input and Filter files
idx=~contains(files,inputName);
files=files(idx);
idx=~contains(files,filterName);
files=files(idx);
%loop through all sample files
for i=1:length(files)
    infile = char(files(i));
    sampleName = infile(strfind(infile,'scarf_')+6:strfind(infile,'.bowtiemap')-1);
    Artist_run_v2(genome,inputName,sampleName,filterName);
end


