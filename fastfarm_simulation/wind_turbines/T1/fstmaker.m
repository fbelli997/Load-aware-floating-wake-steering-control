% fstmaker.m
% Create 8001.fstf ... 8400.fstf by editing line 40 in 5454.fst

srcFile = '15201.fst';
txt = fileread(srcFile);

% Detect original newline sequence to preserve it
newlinePattern = '\r\n|\n|\r';
newlineChar = regexp(txt, newlinePattern, 'match', 'once');
if isempty(newlineChar), newlineChar = sprintf('\n'); end

% Split into lines (cell array of char), preserving empties
lines = regexp(txt, newlinePattern, 'split');
hadTrailingNewline = ~isempty(regexp(txt, [newlinePattern '$'], 'once'));

if numel(lines) < 40
    error('File has only %d lines; cannot modify line 40.', numel(lines));
end

line40 = lines{40};

% Patterns to identify the value for ServoFile on line 40 (try in order):
% 1) Quoted value before the ServoFile label
pat1 = '^(\s*)"[^"]*"(.*\bServoFile\b.*)$';
% 2) Unquoted value before the ServoFile label
pat2 = '^(\s*)(\S+)(\s+\bServoFile\b.*)$';
% 3) First quoted string anywhere on the line (fallback)
pat3 = '"[^"]*"';

for k = 15202:15212
    newQuoted = sprintf('"turbine_controllers/servodyn_%04d.dat"', k);
    newLine = line40; matched = false;

    % Try quoted + ServoFile
    tok = regexp(newLine, pat1, 'tokens', 'once');
    if ~isempty(tok)
        newLine = [tok{1} newQuoted tok{2}];
        matched = true;
    end

    % Try unquoted + ServoFile
    if ~matched
        tok = regexp(newLine, pat2, 'tokens', 'once');
        if ~isempty(tok)
            newLine = [tok{1} newQuoted tok{3}];
            matched = true;
        end
    end

    % Fallback: replace first quoted token anywhere on the line
    if ~matched && ~isempty(regexp(newLine, pat3, 'once'))
        newLine = regexprep(newLine, pat3, newQuoted, 'once');
        matched = true;
    end

    if ~matched
        warning('Could not identify ServoFile value on line 40. Writing unchanged for %d.', k);
    end

    % Assemble output with preserved newlines/spaces
    outLines = lines;
    outLines{40} = newLine;

    outTxt = outLines{1};
    for i = 2:numel(outLines)
        outTxt = [outTxt, newlineChar, outLines{i}]; %#ok<AGROW>
    end
    if hadTrailingNewline
        outTxt = [outTxt, newlineChar]; %#ok<AGROW>
    end

    outName = sprintf('%04d.fst', k);
    fid = fopen(outName, 'w');
    if fid == -1, error('Could not open %s for writing.', outName); end
    fwrite(fid, outTxt, 'char');
    fclose(fid);
end

disp('Done: generated 9001.fstf through 9200.fstf.');
