function pairs = find_opposite_saccades(ms, varargin)
% FIND_OPPOSITE_SACCADES  Find pairs of saccades with similar amplitude,
% opposite direction, and close in time.
%
%   pairs = FIND_OPPOSITE_SACCADES(ms)
%   pairs = FIND_OPPOSITE_SACCADES(ms, 'TolAmp',0.1, 'TolDir',20, 'MaxGap',50)
%
% Inputs:
%   ms         - N-by-8 array with columns:
%                1 startIdx, 2 endIdx, 3 peakVel, 4 amplitude,
%                5 direction (deg), 6 horizAmp, 7 vertAmp, 8 displacement
%
% Optional name/value pairs:
%   'TolAmp'   - relative amplitude tolerance (default 0.10 for 10%)
%                (matching criterion: abs(a1-a2) <= TolAmp * mean([a1,a2]))
%   'TolDir'   - angular tolerance in degrees around 180 (default 20)
%                (direction difference criterion: |wrapTo180(d1-d2)| within
%                 (180 +/- TolDir))
%   'MaxGap'   - maximum allowed gap in samples between end of first and
%                start of second (default 50). Only pairs where start_j > end_i
%                AND (start_j - end_i) <= MaxGap are considered.
%   'ReturnTable' - true/false (default false). If true, returns table; otherwise matrix.
%
% Output:
%   pairs - table with columns:
%           idx1, idx2, start1, end1, start2, end2, amp1, amp2, dirDiffDeg,
%           ampRelDiff, lagSamples
%           (If ReturnTable==false, returns numeric matrix with same columns.)
%
% Example:
%   pairs = find_opposite_saccades(ms, 'TolAmp',0.10, 'TolDir',20, 'MaxGap',100);

% Parse inputs
p = inputParser;
p.addRequired('ms', @(x) isnumeric(x) && size(x,2) >= 8);
p.addParameter('TolAmp', 0.10, @(x) isnumeric(x) && x>=0);
p.addParameter('TolDir', 20, @(x) isnumeric(x) && x>=0 && x<=180);
p.addParameter('MaxGap', 50, @(x) isnumeric(x) && x>=0);
p.addParameter('ReturnTable', false, @(x) islogical(x) || isnumeric(x));
p.parse(ms, varargin{:});
tolAmp = p.Results.TolAmp;
tolDir = p.Results.TolDir;
maxGap = p.Results.MaxGap;
returnTable = logical(p.Results.ReturnTable);

% Extract columns
startIdx = ms(:,1);
endIdx   = ms(:,2);
amp      = ms(:,4);
dirDeg   = ms(:,5);

N = size(ms,1);
results = zeros(0,11); % will collect rows: idx1 idx2 start1 end1 start2 end2 amp1 amp2 dirDiff ampRelDiff lag

% For efficiency: create index list sorted by start time
[~, order] = sort(startIdx);
startIdxS = startIdx(order);
endIdxS   = endIdx(order);
ampS      = amp(order);
dirS      = dirDeg(order);

% Build an index that maps sorted order back to original indices
origIdx = (1:N)';
origIdxS = origIdx(order);

% For each saccade i, find candidate j where start_j > end_i and gap <= maxGap
% because of sorting by start, we can use a moving pointer for j_min for efficiency.
j_start_ptr = 1;
for ii = 1:N
    i_orig = origIdxS(ii);
    i_end = endIdxS(ii);
    i_amp = ampS(ii);
    i_dir = dirS(ii);

    % advance j_start_ptr until start_j > end_i (we want j strictly after i in time)
    while j_start_ptr <= N && startIdxS(j_start_ptr) <= i_end
        j_start_ptr = j_start_ptr + 1;
    end
    if j_start_ptr > N
        continue
    end

    % find last candidate whose start <= i_end + maxGap
    % use binary-style scan (vectorized search)
    j_last = find(startIdxS <= i_end + maxGap, 1, 'last');
    if isempty(j_last) || j_last < j_start_ptr
        continue
    end

    % candidate range
    candIdxs = j_start_ptr:j_last;

    % compute amplitude relative difference: abs(a1-a2) / mean([a1,a2])
    a2 = ampS(candIdxs);
    ampRelDiff = abs(i_amp - a2) ./ ((i_amp + a2)/2 + eps); % avoid divide-by-zero

    ampMask = ampRelDiff <= tolAmp;

    % compute angular difference to 180 deg
    d2 = dirS(candIdxs);
    % compute signed angular difference in degrees in range [-180,180]
    angDiff = wrapTo180(i_dir - d2);       % negative means d2 > i_dir
    angDiffAbs = abs(angDiff);
    % For "opposite", we want angDiffAbs near 180: equivalent to |angDiffAbs - 180| <= tolDir
    angFrom180 = abs(angDiffAbs - 180);
    dirMask = angFrom180 <= tolDir;

    % combined mask
    valid = find(ampMask & dirMask);
    if isempty(valid)
        continue
    end

    % add rows for each valid candidate
    for k = 1:numel(valid)
        jj = candIdxs(valid(k));
        j_orig = origIdxS(jj);
        lag = startIdxS(jj) - i_end; % >=1, <= maxGap
        row = [i_orig, j_orig, startIdxS(ii), endIdxS(ii), startIdxS(jj), endIdxS(jj), ...
               i_amp, ampS(jj), angFrom180(valid(k)), ampRelDiff(valid(k)), lag];
        results(end+1, :) = row; %#ok<AGROW>
    end
end

% Build output
if isempty(results)
    if returnTable
        pairs = table([],[],[],[],[],[],[],[],[],[],[],'VariableNames', ...
            {'idx1','idx2','start1','end1','start2','end2','amp1','amp2','dirFrom180Deg','ampRelDiff','lagSamples'});
    else
        pairs = zeros(0,11);
    end
    return
end

if returnTable
    pairs = table(results(:,1), results(:,2), results(:,3), results(:,4), results(:,5), results(:,6), ...
                  results(:,7), results(:,8), results(:,9), results(:,10), results(:,11), ...
                  'VariableNames', {'idx1','idx2','start1','end1','start2','end2','amp1','amp2','dirFrom180Deg','ampRelDiff','lagSamples'});
else
    pairs = results;
end

end

% --- helper: wrapTo180 (keeps angles in [-180,180])
function a = wrapTo180(x)
    a = mod(x + 180, 360) - 180;
end
