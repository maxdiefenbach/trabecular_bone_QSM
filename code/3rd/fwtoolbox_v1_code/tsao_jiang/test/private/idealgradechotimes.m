%
% [echotimes] = idealgradechotimes(MinEchoTimeInMs, [fieldstrengths]),
%
% IDEAL water fat suppression - gradient echo

% 
% 2006.10.19 Jeffrey Tsao
% 2007.08.24 Jeffrey Tsao
%  - corrected based on Pineda et al. MRM. 54: 625-635. 2005
%    Mid-echo should be 90-phase offset
% 2011.09.16 Jeffrey Tsao - Added waterfatppmdiff
function [echotimes] = idealgradechotimes(MinEchoTimeInMs, fieldstrengths, waterfatppmdiff),
if nargin<1, help(mfilename); end
if nargin<2, fieldstrengths=[]; end
if nargin<3, waterfatppmdiff=[]; end

if isempty(fieldstrengths), fieldstrengths = 7; end % Tesla
if isempty(waterfatppmdiff), waterfatppmdiff = 3.4; end % Water fat difference

GyromagneticRatio = 42.58;         % MHz/T

fprintf('Assumed Larmor frequency:');
fprintf(' %f',fieldstrengths*GyromagneticRatio);
fprintf('\n');
LarmorFrequency = fieldstrengths*GyromagneticRatio;


HzDif = (LarmorFrequency) * waterfatppmdiff;   % Hz/ppm*ppm
MsForOneCycle = 1000./HzDif;
fprintf('Time for one cycle of water-fat difference: %fms\n',MsForOneCycle)
NumTwelfthCycles = MinEchoTimeInMs(1)./(MsForOneCycle/12);   % Min Number of 1/12 cycles
NumTwelfthCycles = ceil(NumTwelfthCycles); 

% Round to next 5 + 6*n with n being an integer
NumTwelfthCycles = ceil((NumTwelfthCycles-5)/6)*6+5;

echotimes = [NumTwelfthCycles;NumTwelfthCycles+4;NumTwelfthCycles+8] * diag(MsForOneCycle/12);

for n=[1:length(LarmorFrequency)],
  fprintf('IDEAL echo times at %.1fT =\n',fieldstrengths(n));
  for m=[1:size(echotimes,1)],
    fprintf(' %.3f (%+.0f/12 cycles)\n',echotimes(m,n),mod(echotimes(m,n)/MsForOneCycle(n)*12 +6,12)-6);
  end
  fprintf('\n');
end; clear n;
