saveDir = "C:/Users/nick/Projects/DATA_STRUCTS";
dataDir = "C:/Users/nick/Projects/RAW_DATA";
codeDir = "C:/Users/nick/Projects/decodingProject";
addpath(genpath(codeDir))

group = fmr1CircTrack_x_makeStruct(dataDir, saveDir);
replayEvents = nick_makeReplayStruct(group, saveDir);
runEvents = nick_makeRunStruct(group, saveDir);
nick_analyzeSync(replayEvents, saveDir);
nick_analyzeWP(group, replayEvents, saveDir);
nick_analyzeRun(runEvents, saveDir);