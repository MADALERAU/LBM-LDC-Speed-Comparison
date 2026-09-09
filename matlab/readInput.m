function [numX,numY,maxItt, initialDensity, reynoldsNumber,tau,maxError] = readInput(filepath)
    lines = readlines(filepath);
    processedLines = cell(6,1);

    for i = 1:length(lines)
        str = lines(i);
        if ~strcmp(str, "")
            splitLines = split(str, '[');
            processedLines{i} = erase(splitLines(2), ']');
        end
    end
    
    numX = str2double(processedLines{1});
    numY = str2double(processedLines{1});
    maxItt = str2double(processedLines{2});
    initialDensity = str2double(processedLines{3});
    reynoldsNumber = str2double(processedLines{4});
    tau = str2double(processedLines{5});
    maxError = str2double(processedLines{6});
end