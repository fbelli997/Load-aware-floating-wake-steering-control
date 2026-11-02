function generate_fstf_files(baseFile, outDir, onlyRunID)
% generate_fstf_files Read base d0.fstf and produce .fstf variants per rules.
% Usage:
%   generate_fstf_files('d0.fstf','fstf_outputs')                % all
%   generate_fstf_files('d0.fstf','fstf_outputs','8045')         % just 8045
%
% If it "runs without error but writes nothing", this version will print
% exactly what it's doing so you can see where it stops.

wind_speed = '9';
if nargin < 1 || isempty(baseFile), baseFile = [wind_speed,'.fstf']; end
if nargin < 2 || isempty(outDir),   outDir   = 'fstf_outputs'; end

fprintf('Base: %s\n', baseFile);
if ~isfile(baseFile)
    error('Base file "%s" not found in %s', baseFile, pwd);
end
if ~exist(outDir,'dir')
    fprintf('Creating output folder: %s\n', outDir);
    mkdir(outDir);
end

% ------- SPEC LIST (from your message) -------
specLines = [
"8010	w7_	d0	nogpr_	gre_	seed1"
"8011	w7_	d0	nogpr_	wst_	seed1"
"8012	w7_	d0	nogpr_	ind_	seed1"
"8013	w7_	d0	nogpr_	hyb_	seed1"
"8014	w8_	d0	nogpr_	gre_	seed1"
"8015	w8_	d0	nogpr_	wst_	seed1"
"8016	w8_	d0	nogpr_	ind_	seed1"
"8017	w8_	d0	nogpr_	hyb_	seed1"
"8018	w9_	d0	nogpr_	gre_	seed1"
"8019	w9_	d0	nogpr_	wst_	seed1"
"8020	w9_	d0	nogpr_	ind_	seed1"
"8021	w9_	d0	nogpr_	hyb_	seed1"
"8022	w7_	dn5	nogpr_	gre_	seed1"
"8023	w7_	dn5	nogpr_	wst_	seed1"
"8024	w7_	dn5	nogpr_	ind_	seed1"
"8025	w7_	dn5	nogpr_	hyb_	seed1"
"8026	w8_	dn5	nogpr_	gre_	seed1"
"8027	w8_	dn5	nogpr_	wst_	seed1"
"8028	w8_	dn5	nogpr_	ind_	seed1"
"8029	w8_	dn5	nogpr_	hyb_	seed1"
"8030	w9_	dn5	nogpr_	gre_	seed1"
"8031	w9_	dn5	nogpr_	wst_	seed1"
"8032	w9_	dn5	nogpr_	ind_	seed1"
"8033	w9_	dn5	nogpr_	hyb_	seed1"
"8034	w7_	dn10	nogpr_	gre_	seed1"
"8035	w7_	dn10	nogpr_	wst_	seed1"
"8036	w7_	dn10	nogpr_	ind_	seed1"
"8037	w7_	dn10	nogpr_	hyb_	seed1"
"8038	w8_	dn10	nogpr_	gre_	seed1"
"8039	w8_	dn10	nogpr_	wst_	seed1"
"8040	w8_	dn10	nogpr_	ind_	seed1"
"8041	w8_	dn10	nogpr_	hyb_	seed1"
"8042	w9_	dn10	nogpr_	gre_	seed1"
"8043	w9_	dn10	nogpr_	wst_	seed1"
"8044	w9_	dn10	nogpr_	ind_	seed1"
"8045	w9_	dn10	nogpr_	hyb_	seed1"
"8046	w7_	dn5	gpr	gre_	seed1"
"8047	w7_	dn5	gpr	wst_	seed1"
"8048	w7_	dn5	gpr	ind_	seed1"
"8049	w7_	dn5	gpr	hyb_	seed1"
"8050	w8_	dn5	gpr	gre_	seed1"
"8051	w8_	dn5	gpr	wst_	seed1"
"8052	w8_	dn5	gpr	ind_	seed1"
"8053	w8_	dn5	gpr	hyb_	seed1"
"8054	w9_	dn5	gpr	gre_	seed1"
"8055	w9_	dn5	gpr	wst_	seed1"
"8056	w9_	dn5	gpr	ind_	seed1"
"8057	w9_	dn5	gpr	hyb_	seed1"
];

% Optional filter: just one RUNID
if nargin >= 3 && ~isempty(onlyRunID)
    mask = false(size(specLines));
    for i = 1:numel(specLines)
        if startsWith(specLines(i), string(onlyRunID))
            mask(i) = true;
        end
    end
    specLines = specLines(mask);
    if isempty(specLines)
        error('No spec line found starting with RUNID "%s".', onlyRunID);
    end
    fprintf('Filtering to RUNID=%s (count=%d)\n', onlyRunID, numel(specLines));
else
    fprintf('Preparing to generate %d files.\n', numel(specLines));
end

% ------- CONSTANT TEXT BLOCKS (exact spacing preserved) -------
if wind_speed == '7'
    blk_dn0 = [
"300      -210     0.000    ""wind_turbines/T1/8001.fst""     150      -360     3       2.8       5       5"
"1500     -70      0.000    ""wind_turbines/T2/8001.fst""     1350     -220     3       2.8       5       5"
"2700      70      0.000    ""wind_turbines/T3/8001.fst""     2550      -80      3      2.8      5       5"
"3900      210     0.000    ""wind_turbines/T4/8001.fst""     3750      60     3        2.8    5       5"
];
    blk_dn5 = [
"300      -50      0.000    ""wind_turbines/T1/8001.fst""     150      -200     3       2.8       5       5"
"1510     -20      0.000    ""wind_turbines/T2/8001.fst""     1360     -170     3       2.8       5       5"
"2715      20      0.000    ""wind_turbines/T3/8001.fst""     2565      -130    3       2.8     5       5"
"3925      50      0.000    ""wind_turbines/T4/8001.fst""     3775      -100    3       2.8       5       5"
];
    blk_dn10 = [
"300       	105       0.000    ""wind_turbines/T1/8001.fst""     150       -45      3       2.8       5       5"
"1505     	35        0.000    ""wind_turbines/T2/8001.fst""     1355      -115     3       2.8       5       5"
"2710     	-35       0.000    ""wind_turbines/T3/8001.fst""     2560      -185     3       2.8       5       5"
"3920      	-105      0.000    ""wind_turbines/T4/8001.fst""     3770      -255     3       2.8       5       5"
];

elseif wind_speed == '8'
        blk_dn0 = [
"300      -210     0.000    ""wind_turbines/T1/8001.fst""     150      -360     3       4       5       5"
"1500     -70      0.000    ""wind_turbines/T2/8001.fst""     1350     -220     3       4       5       5"
"2700      70      0.000    ""wind_turbines/T3/8001.fst""     2550      -80      3      4      5       5"
"3900      210     0.000    ""wind_turbines/T4/8001.fst""     3750      60     3       4    5       5"
];
    blk_dn5 = [
"300      -50      0.000    ""wind_turbines/T1/8001.fst""     150      -200     3       4       5       5"
"1510     -20      0.000    ""wind_turbines/T2/8001.fst""     1360     -170     3       4       5       5"
"2715      20      0.000    ""wind_turbines/T3/8001.fst""     2565      -130    3       4     5       5"
"3925      50      0.000    ""wind_turbines/T4/8001.fst""     3775      -100    3       4       5       5"
];
    blk_dn10 = [
"300       	105       0.000    ""wind_turbines/T1/8001.fst""     150       -45      3       4       5       5"
"1505     	35        0.000    ""wind_turbines/T2/8001.fst""     1355      -115     3       4       5       5"
"2710     	-35       0.000    ""wind_turbines/T3/8001.fst""     2560      -185     3       4       5       5"
"3920      	-105      0.000    ""wind_turbines/T4/8001.fst""     3770      -255     3       4       5       5"
];
elseif wind_speed == '9'
        blk_dn0 = [
"300      -210     0.000    ""wind_turbines/T1/8001.fst""     150      -360     3       4.5       5       5"
"1500     -70      0.000    ""wind_turbines/T2/8001.fst""     1350     -220     3       4.5       5       5"
"2700      70      0.000    ""wind_turbines/T3/8001.fst""     2550      -80      3      4.5      5       5"
"3900      210     0.000    ""wind_turbines/T4/8001.fst""     3750      60     3       4.5    5       5"
];
    blk_dn5 = [
"300      -50      0.000    ""wind_turbines/T1/8001.fst""     150      -200     3       4.5       5       5"
"1510     -20      0.000    ""wind_turbines/T2/8001.fst""     1360     -170     3       4.5       5       5"
"2715      20      0.000    ""wind_turbines/T3/8001.fst""     2565      -130    3       4.5     5       5"
"3925      50      0.000    ""wind_turbines/T4/8001.fst""     3775      -100    3       4.5       5       5"
];
    blk_dn10 = [
"300       	105       0.000    ""wind_turbines/T1/8001.fst""     150       -45      3       4.5       5       5"
"1505     	35        0.000    ""wind_turbines/T2/8001.fst""     1355      -115     3       4.5       5       5"
"2710     	-35       0.000    ""wind_turbines/T3/8001.fst""     2560      -185     3       4.5       5       5"
"3920      	-105      0.000    ""wind_turbines/T4/8001.fst""     3770      -255     3       4.5       5       5"
];
end


% ------- INFLOW MAP (line 34) -------
inflowMap = containers.Map( ...
    {'w7_','w8_','w9_'}, ...
    { ...
    '"inflow/seed_wind7.dat"  InflowFile         - Name of file containing InflowWind module input parameters (quoted string)', ...
    '"inflow/seed_wind8.dat"  InflowFile         - Name of file containing InflowWind module input parameters (quoted string)', ...
    '"inflow/seed_wind9.dat"  InflowFile         - Name of file containing InflowWind module input parameters (quoted string)' ...
    } ...
);

% ------- READ BASE FILE (robust split on \r\n or \n) -------
raw = fileread(baseFile);
baseLines = regexp(raw, '\r\n|\n', 'split');  % cellstr
if isempty(baseLines) || (numel(baseLines)==1 && baseLines{1}=="")
    error('Base file appears empty.');
end
% remove a possible trailing empty last line
if baseLines{end} == "", baseLines(end) = []; end
fprintf('Base file has %d lines.\n', numel(baseLines));
if numel(baseLines) < 42
    error('Base file must have at least 42 lines (has %d).', numel(baseLines));
end

nWritten = 0;

% ------- PROCESS EACH SPEC -------
for i = 1:numel(specLines)
    cols = split(strtrim(specLines(i)));
    if numel(cols) < 6
        warning('Skipping malformed spec: "%s"', specLines(i));
        continue;
    end
    runID    = cols{1};
    wToken   = cols{2};  % w7_/w8_/w9_
    dnToken  = cols{3};  % d0/dn5/dn10
    gprToken = cols{4};  % gpr/nogpr_
    dirToken = cols{5};  % gre_/wst_/ind_/hyb_
    seed     = cols{6};  % seed1

    lines = baseLines;   % cellstr

    % --- line 34 (1-based) ---
    if ~isKey(inflowMap, wToken)
        error('Unknown wind token: %s (spec: %s)', wToken, specLines(i));
    end
    lines{34} = inflowMap(wToken);

    % --- lines 39..42 (1-based) ---
    switch lower(dnToken)
        case {'d0','dn0'}
            blk = blk_dn0;
        case 'dn5'
            blk = blk_dn5;
        case 'dn10'
            blk = blk_dn10;
        otherwise
            error('Unknown dn token: %s (spec: %s)', dnToken, specLines(i));
    end
    blk = replace(blk, "8001", runID);  % string array
    % convert to cellstr for uniform writing
    lines(39:42) = cellstr(blk);

    % --- filename ---
    outName = makeFilename(runID, wToken, dnToken, gprToken, dirToken, seed);
    outPath = fullfile(outDir, outName);

    % --- write file ---
    fid = fopen(outPath,'w');
    if fid < 0
        error('Cannot write file: %s', outPath);
    end
    for k = 1:numel(lines)
        fprintf(fid, '%s\n', char(lines{k}));
    end
    fclose(fid);
    nWritten = nWritten + 1;
    fprintf('Wrote %s\n', outPath);
end

fprintf('Done. %d files written to: %s\n', nWritten, outDir);

end % function

function s = stripUnderscore(s)
if endsWith(s, "_"), s = extractBefore(s, strlength(s)); end
end

function fname = makeFilename(runID,wToken,dnToken,gprToken,dirToken,seed)
parts = { ...
    runID, ...
    %stripUnderscore(wToken), ...
    %dnToken, ...
    %stripUnderscore(gprToken), ...
    %stripUnderscore(dirToken), ...
    %seed ...
    };
fname = strjoin(parts, "_") + ".fstf";
fname = char(fname); % ensure char for fullfile/fopen on older MATLAB
end
