%%
% failed = VerifyEchoTimes(TEs_Seconds, FieldStrength_Tesla, WaterFat_ppmdif, verbose)
%
% Verify if given echo times are at optimal settings for water and fat to
% be phased according to Pineda et al. Magn Reson Med. 2005 Sep;54(3):625-35.
%
% Output: failed
%  WARNING
%  -1 - Phases of 1st and 3rd TE are swapped.
%  -2 - Phase of fat signal at 2nd TE is not optimal.
%  ERROR
%  1 - number of TE is not 3
%  2 - Phase of fat signal at 1st TE is not 60 degree from that at 2nd TE.
%  3 - Phase of fat signal at 3rd TE is not 60 degree from that at 2nd TE.
%  4 - Fat signal at 1st and 3rd TE is at the same phase.
%  

% Jeffrey Tsao, 2011
function failed = VerifyEchoTimes(TEs_Seconds, FieldStrength_Tesla, WaterFat_ppmdif, verbose)
    failed = 0;
    if nargin<4, verbose=[]; end
    if isempty(verbose), verbose=1; end
    if isempty(WaterFat_ppmdif) WaterFat_ppmdif = 3.4; end          % Water fat chemical shift difference in ppm

    if numel(TEs_Seconds)~=3,
        if verbose, fprintf('Number of TEs should be 3.'); end
        failed = 1; return;
    end
    GyromagneticRatio = 42.576;         % MHz/T
    LarmorFrequency = FieldStrength_Tesla*GyromagneticRatio;
    HzDif = (LarmorFrequency) * WaterFat_ppmdif;   % Hz/ppm*ppm
    SecondForOneCycle = 1.0/HzDif;
    NumTwelfthCycles = mod(round(TEs_Seconds(:) / SecondForOneCycle * 12.0),12);

    % First or third echo should be at 4/12 of a cycle away from middle echo
    dif_1 = mod(NumTwelfthCycles(1)-NumTwelfthCycles(2),12);
    dif_3 = mod(NumTwelfthCycles(3)-NumTwelfthCycles(2),12);
    if dif_1~=4 && dif_1~=8,
        if verbose,
            fprintf('Echos at %.0f/12 %.0f/12 %.0f/12 of circle\nThird echo is incorrect.\n',NumTwelfthCycles(1:3));
        end
        failed = 2; return;
    end
    if dif_3~=4 && dif_3~=8,
        if verbose,
            fprintf('Echos at %.0f/12 %.0f/12 %.0f/12 of circle\nThird echo is incorrect.\n',NumTwelfthCycles(1:3));
        end
        failed = 3; return;
    end

    % First and third echo should not be the same
    if dif_1==dif_3,
        if verbose,
            fprintf('Echos at %.0f/12 %.0f/12 %.0f/12 of circle\nFirst and third echo are the same.\n',NumTwelfthCycles(1:3));
        end
        failed = 4; return;
    end

    % Middle echo should be at 3/12 or 9/12 of a cycle
    if NumTwelfthCycles(2)~=3 && NumTwelfthCycles(2)~=9,
        if verbose,
            fprintf('WARNING: Middle echo is not optimal.\n         Water/fat phase differences at %.0f/12 %.0f/12 %.0f/12 of unit circle.\n',NumTwelfthCycles(1:3));
        end
        failed = -2;
    end

    % If first and third echo are swapped
    if dif_1==4 && dif_3==8, 
        fprintf('WARNING: 1st and third echoes are swapped.\n');
        failed = -1;   % This needs to be after than failed = -2, because it is more serious (i.e. need to be corrected).
    end

end
