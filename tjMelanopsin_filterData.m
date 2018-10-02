function dataTraceFiltered = tjMelanopsin_interpolateData(dataTraceInterpolated)
dataTraceFiltered = sgolayfilt(dataTraceInterpolated, 11, 21);