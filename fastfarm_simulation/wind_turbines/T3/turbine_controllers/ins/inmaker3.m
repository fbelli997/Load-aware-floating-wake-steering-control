% zmqmaker.m
% Create I8001.IN ... I8400.IN by editing line 187 in IEA15MW_ZeroMQ.IN

srcFile = 'I15001.IN';
txt = fileread(srcFile);

% Detect original newline sequence to preserve it
newlinePattern = '\r\n|\n|\r';
newlineChar = regexp(txt, newlinePattern, 'match', 'once');
if isempty(newlineChar), newlineChar = sprintf('\n'); end

% Split into lines (cell array of char), preserving empties
lines = regexp(txt, newlinePattern, 'split');
hadTrailingNewline = ~isempty(regexp(txt, [newlinePattern '$'], 'once'));

if numel(lines) < 187
    error('File has only %d lines; cannot modify line 187.', numel(lines));
end

line187 = lines{187};

% Patterns (try in order):
% 1) Quoted tcp address before the ZMQ_CommAddress label
pat1 = '^(\s*)"tcp://[^"]*"(.*\bZMQ_CommAddress\b.*)$';
% 2) Any quoted token before the ZMQ_CommAddress label (fallback)
pat2 = '^(\s*)"[^"]*"(.*\bZMQ_CommAddress\b.*)$';
% 3) Replace the first quoted tcp address anywhere on the line (last resort)
pat3 = '"tcp://[^"]*"';

for k = 15201:15212
    newQuoted = sprintf('"tcp://localhost:%d"', k);
    newLine = line187; matched = false;

    % Try tcp-specific + label
    tok = regexp(newLine, pat1, 'tokens', 'once');
    if ~isempty(tok)
        newLine = [tok{1} newQuoted tok{2}];
        matched = true;
    end

    % Try generic quoted + label
    if ~matched
        tok = regexp(newLine, pat2, 'tokens', 'once');
        if ~isempty(tok)
            newLine = [tok{1} newQuoted tok{2}];
        matched = true;
        end
    end

    % Replace first quoted tcp address anywhere
    if ~matched && ~isempty(regexp(newLine, pat3, 'once'))
        newLine = regexprep(newLine, pat3, newQuoted, 'once');
        matched = true;
    end

    if ~matched
        warning('Could not identify the address on line 187. Writing unchanged for %d.', k);
    end

    % Write file
    outLines = lines;
    outLines{187} = newLine;

    % Reassemble text with original newline style
    outTxt = outLines{1};
    for i = 2:numel(outLines), outTxt = [outTxt, newlineChar, outLines{i}]; %#ok<AGROW>
    end
    if hadTrailingNewline, outTxt = [outTxt, newlineChar]; end %#ok<AGROW>

    outName = sprintf('I%04d.IN', k);
    fid = fopen(outName, 'w');
    if fid == -1, error('Could not open %s for writing.', outName); end
    fwrite(fid, outTxt, 'char');
    fclose(fid);
end

disp('Done: generated I8001.IN through I8400.IN.');
