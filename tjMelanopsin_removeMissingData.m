function dataTraceRaw = tjMelanopsin_processData(dataTraceRaw0)

lengthVector = length(dataTraceRaw0);
padSz = 5;
dataTraceRaw = dataTraceRaw0;
rmIdxMOD = [find(dataTraceRaw0 == 0) ...
    find(dataTraceRaw0 == 0)-padSz ...
    find(dataTraceRaw0 == 0)+padSz];
rmIdxMOD(rmIdxMOD < 1) = [];
dataTraceRaw(rmIdxMOD) = NaN;

% Truncate with original length
dataTraceRaw = dataTraceRaw(1:lengthVector);