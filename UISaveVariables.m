%% Script that Saves Files on Request %%

VAR = evalin('base','VAR');
LName    = [sprintf('Name ; '     ),sprintf('%s ; ',VAR.Name{:})  ]; LName   = LName(1:end-3);
LHandle  = [sprintf('Handle ; '   ),sprintf('%d ; ',VAR.Handle(:))]; LHandle = LHandle(1:end-3);
LMag     = [sprintf('Magnitude ; '),sprintf('%d ; ',VAR.Mag(:))   ]; LMag    = LMag(1:end-3);
LDat     = [sprintf('Data ; '     ),sprintf('%d ; ',VAR.Dat(:))   ]; LDat    = LDat(1:end-3);
LSim     = [                        sprintf('%s ; ',VAR.Sim{:})   ]; LSim    = LSim(1:end-3);
LPlot    = [                        sprintf('%s ; ',VAR.Plot{:})  ]; LPlot   = LPlot(1:end-3);
LLLP     = [                        sprintf('%s ; ',VAR.LLP{:})   ]; LLLP    = LLLP(1:end-3);
FULLTEXT =  sprintf([LName,'\n',LHandle,'\n',LMag,'\n',LDat,'\n',LSim,'\n',LPlot,'\n',LLLP]);
dlmwrite([PathName,FileName],FULLTEXT,'delimiter','')

