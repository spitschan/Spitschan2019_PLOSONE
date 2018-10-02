function [dataTraceInterpolated, t] = tjMelanopsin_interpolateData(timeTraceRawExpt, dataTraceRawExpt, frequencyHz)
dt = 1/frequencyHz;
t = 0:dt:(timeTraceRawExpt(end)+dt);
dataTraceInterpolated = interp1(timeTraceRawExpt, dataTraceRawExpt, t);