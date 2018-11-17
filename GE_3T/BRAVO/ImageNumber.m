function number = ImageNumber(pass, slice, echo, pfile)
% Image numbering scheme (P = Phase; S = Slice; E = Echo):
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn

    % Need to map the legacy "pass" number to a phase number
    numPassesPerPhase = fix(pfile.passes / pfile.phases);
    phase = fix(pass / numPassesPerPhase);

    slicesPerPhase = pfile.slicesPerPass * numPassesPerPhase * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1) + 1;
end

