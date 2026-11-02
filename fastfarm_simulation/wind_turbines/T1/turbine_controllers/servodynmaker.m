% servodynmaker_v2.m
srcFile = 'IEA-22-280-RWT_ServoDyn.dat';
txt = fileread(srcFile);

% Detect original newline sequence to preserve it
newlinePattern = '\r\n|\n|\r';
newlineChar = regexp(txt, newlinePattern, 'match', 'once');
if isempty(newlineChar), newlineChar = sprintf('\n'); end

% Split into lines (cell array of char), preserving empties
lines = regexp(txt, newlinePattern, 'split');
hadTrailingNewline = ~isempty(regexp(txt, [newlinePattern '$'], 'once'));

if numel(lines) < 78
    error('File has only %d lines; cannot modify line 78.', numel(lines));
end

line78 = lines{78};

% We’ll try these patterns in order:
% 1) Quoted value before DLL_InFile label
pat1 = '^(\s*)"([^"]*)"(\s+DLL_InFile\b.*)$';
% 2) Unquoted value before DLL_InFile label
pat2 = '^(\s*)(\S+)(\s+DLL_InFile\b.*)$';
% 3) Any first quoted string on the line (fallback)
pat3 = '^(\s*)"[^"]*"(.*)$';

for k = 15201:15212
    newQuoted = sprintf('"ins/I%04d.IN"', k);
    newLine = line78; matched = false;

    % Try quoted + DLL_InFile
    tok = regexp(newLine, pat1, 'tokens', 'once');
    if ~isempty(tok)
        newLine = [tok{1} newQuoted tok{3}];
        matched = true;
    end

    % Try unquoted + DLL_InFile
    if ~matched
        tok = regexp(newLine, pat2, 'tokens', 'once');
        if ~isempty(tok)
            newLine = [tok{1} newQuoted tok{3}];
            matched = true;
        end
    end

    % Fallback: replace first quoted token anywhere on the line
    if ~matched
        if ~isempty(regexp(newLine, '"[^"]*"', 'once'))
            newLine = regexprep(newLine, '"[^"]*"', newQuoted, 'once');
            matched = true;
        end
    end

    if ~matched
        warning('Could not identify the value to replace on line 78. Writing files unchanged for k=%d.', k);
    end

    % Write file
    outLines = lines;
    outLines{78} = newLine;

    % Reassemble with original newline
    outTxt = outLines{1};
    for i = 2:numel(outLines), outTxt = [outTxt, newlineChar, outLines{i}]; %#ok<AGROW>
    end
    if hadTrailingNewline, outTxt = [outTxt, newlineChar]; end %#ok<AGROW>

    outName = sprintf('servodyn_%04d.dat', k);
    fid = fopen(outName, 'w');
    if fid == -1, error('Could not open %s for writing.', outName); end
    fwrite(fid, outTxt, 'char');
    fclose(fid);
end

disp('Done: generated servodyn_8001.dat through servodyn_8400.dat.');
