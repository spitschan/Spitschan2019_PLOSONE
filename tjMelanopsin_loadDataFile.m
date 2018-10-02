function [theData theT theDataIdx] = tjMelanopsin_loadDataFile(dataPath)
Mtmp = readtable(dataPath);
theData = Mtmp.diameter;
theDataIdx = Mtmp.index;
theT = Mtmp.timestamp;