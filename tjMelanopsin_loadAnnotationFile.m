function annotIdx = tjMelanopsin_loadAnnotationFile(annotPath)
Mtmp = readtable(annotPath);
annotIdx = Mtmp.index;