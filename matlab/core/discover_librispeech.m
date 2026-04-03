function corpus = discover_librispeech(root_dirs)
% discover_librispeech - recursively find all .flac files under root_dirs.
%
% usage:
%   corpus = discover_librispeech('/path/to/dev-clean');
%   corpus = discover_librispeech({'/path/to/dev-clean', '/path/to/dev-other'});
%
% returns a struct array with fields:
%   .path     - full path to .flac file
%   .label    - speakerID_chapterID_uttID
%   .speaker  - speaker folder name
%   .chapter  - chapter folder name
%   .corpus   - source corpus tag (e.g. 'dev-clean')

if ischar(root_dirs) || isstring(root_dirs)
    root_dirs = {char(root_dirs)};
end

corpus = struct('path',{},'label',{},'speaker',{},'chapter',{},'corpus',{});

for rd = 1:numel(root_dirs)
    this_root = root_dirs{rd};

    parts = strsplit(strrep(this_root, '\', '/'), '/');
    parts = parts(~cellfun(@isempty, parts));
    corpus_tag = parts{end};

    if ~isfolder(this_root)
        warning('discover_librispeech: not found: %s', this_root);
        continue;
    end

    fprintf('Scanning: %s\n', this_root);
    tic;
    all_files = recursive_find(this_root, '.flac');

    if isempty(all_files)
        warning('discover_librispeech: no .flac files under %s', this_root);
        continue;
    end

    n_before = numel(corpus);
    for i = 1:numel(all_files)
        fpath = all_files{i};
        [fdir, fname, ~] = fileparts(fpath);
        label = strrep(fname, '-', '_');

        dir_parts = strsplit(strrep(fdir, '\', '/'), '/');
        dir_parts = dir_parts(~cellfun(@isempty, dir_parts));
        np = numel(dir_parts);
        if np >= 2
            chapter = dir_parts{np};
            speaker = dir_parts{np-1};
        elseif np == 1
            chapter = dir_parts{1}; speaker = 'unknown';
        else
            chapter = 'unknown'; speaker = 'unknown';
        end

        corpus(end+1).path    = fpath;   %#ok<AGROW>
        corpus(end).label     = label;
        corpus(end).speaker   = speaker;
        corpus(end).chapter   = chapter;
        corpus(end).corpus    = corpus_tag;
    end

    n_added    = numel(corpus) - n_before;
    n_speakers = numel(unique({corpus(n_before+1:end).speaker}));
    fprintf('  %d clips from %d speakers (%.1fs)\n', n_added, n_speakers, toc);
end

fprintf('Total corpus: %d clips.\n\n', numel(corpus));
end


function files = recursive_find(root, ext)
files = {};
items = dir(root);
for i = 1:numel(items)
    if items(i).name(1) == '.', continue; end
    full = fullfile(root, items(i).name);
    if items(i).isdir
        sub = recursive_find(full, ext);
        files = [files, sub]; %#ok<AGROW>
    elseif endsWith_ci(items(i).name, ext)
        files{end+1} = full; %#ok<AGROW>
    end
end
end


function tf = endsWith_ci(str, ext)
tf = strcmpi(str(max(1,end-numel(ext)+1):end), ext);
end
