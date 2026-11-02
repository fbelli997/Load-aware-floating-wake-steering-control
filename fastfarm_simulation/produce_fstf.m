inputFile = '14020.fstf';
lines = readlines(inputFile);

targetLines = 39:42;
startNum = 14021;
endNum   = 14036;

for num = startNum:endNum
    newLines = lines;

    % Change line 34 if needed
    %newLines(34) = sprintf('"inflow/wind%d.bts"\t\t InflowFile         - Name of file containing InflowWind module input parameters (quoted string)', num);

    % Change turbine .fst file numbers in lines 39:42
    for k = targetLines
        newLines(k) = regexprep(newLines(k), '\d+\.fst', sprintf('%d.fst', num));
    end

    outputFile = sprintf('%d.fstf', num);
    writelines(newLines, outputFile);
end