% e5_iso_budget - wilcoxon signed-rank significance tests.
%
% compares uniform and linear against proposed on distortion and stoi.
% reads per-clip csvs from E1 and E6.
%
% outputs: E5_iso_wilcoxon_distortion.csv, E5_iso_wilcoxon_stoi.csv
%
% run after E1 and E6.
clear; clc; close all;
DATA_DIR = './';

if ~exist('signrank','file')
    error('Requires Statistics and Machine Learning Toolbox.'); end

BUDGET_PCTS = [40,60,80,100];
STRATEGIES  = {'Uniform','Linear','Quadratic'};
REGIONS     = {'Mid-Band(1-4kHz)','High-Band(>4kHz)'};
comp_names  = {'Uniform','Linear'};   % compared against Quadratic

% ── Load E1 distortion data ───────────────────────────────────────────────
fprintf('Loading E1 distortion data...\n');
dist_data = struct();
for b=1:numel(BUDGET_PCTS)
    pct=BUDGET_PCTS(b);
    fname=fullfile(DATA_DIR,sprintf('E1_iso_envelope_correlation_%dpct.csv',pct));
    if ~isfile(fname), error('Missing %s — run E1 first.',fname); end
    T=readtable(fname,'TextType','string','VariableNamingRule','preserve');
    clips=unique(T.Clip,'stable'); N=numel(clips);
    mid_mat=zeros(N,3); high_mat=zeros(N,3);   % 3 strategies
    for c=1:N
        for r=1:2
            mask=strcmp(T.Clip,clips(c))&strcmp(T.Region,REGIONS{r});
            row=T(mask,:); if isempty(row), continue; end
            vals=table2array(row(:,5:7));  % Uniform, Linear, Quadratic
            if r==1, mid_mat(c,:)=vals(1,:); else, high_mat(c,:)=vals(1,:); end
        end
    end
    field=sprintf('b%d',pct);
    dist_data.(field).mid=mid_mat; dist_data.(field).high=high_mat;
    dist_data.(field).N=N;
    fprintf('  %d%%: %d clips\n',pct,N);
end

% ── Load E6 STOI data ─────────────────────────────────────────────────────
fprintf('Loading E6 STOI data...\n');
stoi_data=struct();
for b=1:numel(BUDGET_PCTS)
    pct=BUDGET_PCTS(b);
    fname=fullfile(DATA_DIR,sprintf('E6_iso_stoi_%dpct.csv',pct));
    if ~isfile(fname), error('Missing %s — run E6 first.',fname); end
    T=readtable(fname,'TextType','string','VariableNamingRule','preserve');
    stoi_data.(sprintf('b%d',pct)).scores=table2array(T(:,4:6));  % U, L, Q
    stoi_data.(sprintf('b%d',pct)).N=height(T);
    fprintf('  %d%%: %d clips\n',pct,height(T));
end
fprintf('\n');

wil_dist={}; wil_stoi={};

for b=1:numel(BUDGET_PCTS)
    pct=BUDGET_PCTS(b); field=sprintf('b%d',pct);
    for ri=1:2
        if ri==1, mat=dist_data.(field).mid; else, mat=dist_data.(field).high; end
        quad=mat(:,3);   % column 3 = Quadratic
        for ci=1:2       % compare Uniform (1) and Linear (2) vs Quadratic
            cv=mat(:,ci); ok=~isnan(cv)&~isnan(quad);
            cv=cv(ok); qv=quad(ok);
            [p,h]=signrank(cv,qv,'alpha',0.05);
            lc=log10(cv+eps); lq=log10(qv+eps);
            dl=lc-lq; ps=std([lc;lq]); cohd=mean(dl)/(ps+eps);
            mr=median(cv)/(median(qv)+eps);
            wil_dist{end+1}={pct,REGIONS{ri},comp_names{ci},sum(ok),...
                p,h,cohd,median(cv),iqr(cv),median(qv),iqr(qv),mr};
        end
    end
    mat=stoi_data.(field).scores; quad=mat(:,3);
    for ci=1:2
        cv=mat(:,ci); ok=~isnan(cv)&~isnan(quad); cv=cv(ok); qv=quad(ok);
        [p,h]=signrank(cv,qv,'alpha',0.05);
        ps=std([cv;qv]); cohd=mean(cv-qv)/(ps+eps);
        wil_stoi{end+1}={pct,comp_names{ci},sum(ok),p,h,cohd,...
            median(cv),iqr(cv),median(qv),iqr(qv),median(cv)-median(qv)};
    end
end

fid=fopen('E5_iso_wilcoxon_distortion.csv','w');
fprintf(fid,'BudgetPct,Region,Competitor,N,p_value,Significant,CohensD,Median_Comp,IQR_Comp,Median_Quad,IQR_Quad,MedianRatio\n');
for i=1:numel(wil_dist); r=wil_dist{i};
    fprintf(fid,'%d,%s,%s,%d,%.6e,%d,%.4f,%.8f,%.8f,%.8f,%.8f,%.4f\n',...
            r{1},r{2},r{3},r{4},r{5},r{6},r{7},r{8},r{9},r{10},r{11},r{12}); end
fclose(fid); fprintf('Saved: E5_iso_wilcoxon_distortion.csv\n');

fid=fopen('E5_iso_wilcoxon_stoi.csv','w');
fprintf(fid,'BudgetPct,Competitor,N,p_value,Significant,CohensD,Median_Comp,IQR_Comp,Median_Quad,IQR_Quad,MedianDelta\n');
for i=1:numel(wil_stoi); r=wil_stoi{i};
    fprintf(fid,'%d,%s,%d,%.6e,%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f\n',...
            r{1},r{2},r{3},r{4},r{5},r{6},r{7},r{8},r{9},r{10},r{11}); end
fclose(fid); fprintf('Saved: E5_iso_wilcoxon_stoi.csv\n');

fprintf('\n=== HIGH-BAND DISTORTION SIGNIFICANCE ===\n');
fprintf('%-5s  %-10s  %10s  %4s  %8s  %10s\n','Bud','Competitor','p-value','Sig','CohensD','MedRatio');
fprintf('%s\n',repmat('-',1,55));
for i=1:numel(wil_dist); r=wil_dist{i};
    if strcmp(r{2},'High-Band(>4kHz)')
        fprintf('%-5d  %-10s  %10.3e  %4d  %8.3f  %10.2fx\n',...
                r{1},r{3},r{5},r{6},r{7},r{12}); end; end