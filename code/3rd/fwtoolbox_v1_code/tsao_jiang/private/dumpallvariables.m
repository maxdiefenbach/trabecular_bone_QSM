% dumpallvariables

% April 21, 2008 Jeffrey Tsao
% - Send all variables to console
function dumpallvariables()

VariableNames = evalin('caller','who');
for VariableNum = [1:length(VariableNames)],
  assignin('base',VariableNames{VariableNum},evalin('caller',VariableNames{VariableNum}));
end; clear VariableNum;
clear VariableNames;