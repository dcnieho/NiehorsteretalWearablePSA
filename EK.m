% Part of GazeParser package.
% Copyright (C) 2012-2015 Hiroyuki Sogo.
% Distributed under the terms of the GNU General Public License (GPL).

% Detecting microsaccades with the method of Engbert & Kliegl (2003).
% The class MicroSacc and the binsacc function are port of Engbert and Kliel's
% MATLAB scripts.

% Example
% ---------
% ::
%     microsacc = EK(data(1).L, 250);
%
%
% REFERENCE:
%  Engbert, R., & Kliegl, R. (2003). Microsaccades uncover the orientation
%  of covert attention. Vision Res, 43(9), 1035-1045.
%
% An nx8 array which holds detected microsaccades
%
% ====== ==================================================
% column description
% ====== ==================================================
% 1      index of the starting point of the microsaccade
% 2      index of the termination point of the microsaccade
% 3      peak velocity
% 4      amplitude
% 5      direction
% 6      horizontal amplitude
% 7      vertical amplitude
% 8      displacement
% ====== ==================================================
classdef EK
    properties
        samplingFreq
        velocityType
        vfac
        minSamples
        ms
        radiusX
        radiusY
        v
        vl
        vr
        radiusXl
        radiusYl
        radiusXr
        radiusYr
        msL_mono
        msR_mono
        labels
    end
    
    methods
        function obj = EK(data, samplingFreq, velocityType, vfac, minSamples, min_separation_ms, minAmp, shuffle_velocities, v)
            if nargin < 9
                v = [];
            end
            if nargin < 8 || isempty(shuffle_velocities)
                shuffle_velocities = false;
            end
            if nargin < 7 || isempty(minAmp)
                minAmp = 0;
            end
            if nargin < 6 || isempty(min_separation_ms)
                min_separation_ms = 20;
            end
            if nargin < 5 || isempty(minSamples)
                minSamples = 3;
            end
            if nargin < 4 || isempty(vfac)
                vfac = 6;
            end
            if nargin < 3 || isempty(velocityType)
                velocityType = 'fast';
            end
            
            if ~ismember(velocityType, {'slow', 'fast', 'slow1000', 'slow1001'})
                error('Invalid velocityType.');
            end
            
            minSamples = ceil(minSamples / 1000.0 * samplingFreq);
            min_separation_ms = ceil(min_separation_ms / 1000.0 * samplingFreq);
            
            obj.samplingFreq = samplingFreq;
            obj.velocityType = velocityType;
            obj.vfac = vfac;
            obj.minSamples = minSamples;
            ms = [];
            
            if size(data, 2) == 2
                if isempty(v)
                    v = obj.vecvel(data);
                    if shuffle_velocities
                        v = v(randperm(size(v, 1)), :);
                    end
                    obj.v = v;
                end
                
                [ms, obj.radiusX, obj.radiusY] = obj.microsacc(data, v);
                
                if ~isempty(ms)
                    ms = obj.post_processing(ms, data, v, min_separation_ms);
                end
                
                obj.ms = ms;
            else
                vl = obj.vecvel(data(:, 1:2));
                vr = obj.vecvel(data(:, 3:4));
                obj.vl = vl;
                obj.vr = vr;
                
                [lms, obj.radiusXl, obj.radiusYl] = obj.microsacc(data(:, 1:2), vl);
                [rms, obj.radiusXr, obj.radiusYr] = obj.microsacc(data(:, 3:4), vr);
                
                if ~isempty(lms) && ~isempty(rms)
                    [ms, obj.msL_mono, obj.msR_mono, obj.labels] = binsacc(lms, rms);
                    
                    if ~isempty(ms)
                        x = mean([data(:, 1), data(:, 3)], 2,'omitnan');
                        y = mean([data(:, 2), data(:, 4)], 2,'omitnan');
                        vx = mean([vl(:, 1), vr(:, 1)], 2,'omitnan');
                        vy = mean([vl(:, 2), vr(:, 2)], 2,'omitnan');
                        
                        data = [x, y];
                        v = [vx, vy];
                        obj.v = v;
                        ms = obj.post_processing(ms, data, v, min_separation_ms);
                    end
                end
                
                obj.ms = ms;
            end

            % throw out saccades that are too small
            obj.ms(obj.ms(:,4)<minAmp,:) = [];
        end
        
        function v = vecvel(obj, data)
            N = size(data, 1);
            v = zeros(N, 2);
            if strcmp(obj.velocityType, 'fast')
                v(2:end-1, :) = obj.samplingFreq / 2.0 * (data(3:end, :) - data(1:end-2, :));
            elseif strcmp(obj.velocityType, 'slow')
                v(3:end-2, :) = obj.samplingFreq / 6.0 * (data(5:end, :) + data(4:end-1, :) - data(2:end-3, :) - data(1:end-4, :));
            elseif strcmp(obj.velocityType, 'slow1000')
                v(4:end-3, :) = obj.samplingFreq / 6.0 * (data(7:end, :) - data(1:end-6, :));
            elseif strcmp(obj.velocityType, 'slow1001')
                v(4:end-3, :) = obj.samplingFreq / 12.0 * (data(7:end, :) + data(6:end-1, :) + data(5:end-2, :) - data(3:end-4, :) - data(2:end-5, :) - data(1:end-6, :));
            end
        end
        
        function [sac, radiusX, radiusY] = microsacc(obj, data, vdata)
            msdx = sqrt(median(vdata(:, 1).^2,'omitnan') - median(vdata(:, 1),'omitnan').^2);
            msdy = sqrt(median(vdata(:, 2).^2,'omitnan') - median(vdata(:, 2),'omitnan').^2);
            radiusX = obj.vfac * msdx;
            radiusY = obj.vfac * msdy;
            
            idx = find(((vdata(:, 1) / radiusX).^2 + (vdata(:, 2) / radiusY).^2) > 1);
            
            N = length(idx);
            sac = zeros(0,8);
            if isempty(idx)
                return
            end
            nsac = 0;
            dur = 1;
            a = 1;
            k = 1;
            
            while k < N
                if idx(k+1) - idx(k) == 1
                    dur = dur + 1;
                else
                    if dur >= obj.minSamples
                        b = k;
                        sac = [sac; idx(a), idx(b) + 1]; %#ok<AGROW>
                        nsac = nsac + 1;
                    end
                    a = k + 1;
                    dur = 1;
                end
                k = k + 1;
            end
            
            if dur >= obj.minSamples
                b = k;
                sac = [sac; idx(a), idx(b) + 1];
                nsac = nsac + 1;
            end
            
            if nsac > 0
                sac = [sac, zeros(nsac, 6)];
                vabs = sqrt(vdata(:, 1).^2 + vdata(:, 2).^2);
                for s = 1:size(sac, 1)
                    a = sac(s, 1);
                    b = sac(s, 2);
                    vpeak = max(vabs(a:b));
                    
                    displ = norm(data(a, 1:2) - data(b, 1:2));
                    ampl = max(pdist(data(a:b, 1:2), 'euclidean'));
                    
                    sac(s, 3) = vpeak;
                    sac(s, 4) = ampl;
                    delx = data(b, 1) - data(a, 1);
                    dely = data(b, 2) - data(a, 2);
                    phi = 180 / pi * atan2(dely, delx);
                    sac(s, 5) = phi;
                    sac(s, 6) = delx;
                    sac(s, 7) = dely;
                    sac(s, 8) = displ;
                end
            end
        end
        
        function ms = post_processing(obj, ms, data, vdata, min_separation)
            vabs = sqrt(vdata(:, 1).^2 + vdata(:, 2).^2);
            
            ms_new = ms(1, :);
            for i = 1:size(ms, 1) - 1
                offset1 = ms(i, 2);
                onset2 = ms(i + 1, 1);
                if (onset2 - offset1) > min_separation
                    ms_new = [ms_new; ms(i + 1, :)]; %#ok<AGROW>
                else
                    sac = zeros(1, 8);
                    a = ms_new(end, 1);
                    b = ms(i + 1, 2);
                    
                    sac(1) = a;
                    sac(2) = b;
                    vpeak = max(vabs(a:b));
                    
                    displ = norm(data(a, 1:2) - data(b, 1:2));
                    ampl = max(pdist(data(a:b, 1:2), 'euclidean'));
                    
                    sac(3) = vpeak;
                    sac(4) = ampl;
                    delx = data(b, 1) - data(a, 1);
                    dely = data(b, 2) - data(a, 2);
                    phi = 180 / pi * atan2(dely, delx);
                    sac(5) = phi;
                    sac(6) = delx;
                    sac(7) = dely;
                    sac(8) = displ;
                    ms_new(end, :) = sac;
                end
            end
            
            ms = ms_new;
        end
    end
end


function [sac, L_mono, R_mono, labels] = binsacc(L, R)
    % BINSACC
    % :param L: left eye's microsaccades, detected by MicroSacc.
    % :param R: right eye's microsaccades, detected by MicroSacc.
    
    labels = []; % Encodes hit-miss.

    NL = size(L, 1);
    NR = size(R, 1);

    L_mono = [];
    R_mono = [];
    sac = [];

    % If there are no microsaccades in the left eye
    if NL == 0
        R_mono = R;
        labels = [ones(NR, 1), zeros(NR, 1)];
        sac = [];
        return;
    end

    % If there are no microsaccades in the right eye
    if NR == 0
        L_mono = L;
        labels = [ones(NL, 1), zeros(NL, 1)];
        sac = [];
        return;
    end

    nsac = 0;
    for i = 1:NL
        l1 = L(i, 1); % onset
        l2 = L(i, 2); % offset
        overlap = find((R(:, 2) >= l1) & (R(:, 1) <= l2));
        if ~isempty(overlap)
            r1 = R(overlap(1), 1);
            r2 = R(overlap(1), 2);
            vl = L(i, 3);
            vr = R(overlap(1), 3);
            ampl = L(i, 4);
            ampr = R(overlap(1), 4);
            dxl = L(i, 6);
            dyl = L(i, 7);
            dxr = R(overlap(1), 6);
            dyr = R(overlap(1), 7);
            displ = L(i, 8);
            dispr = R(overlap(1), 8);
            dx = dxl + dxr;
            dy = dyl + dyr;
            phi = 180 / pi * atan2(dy, dx);
            sac = [sac; min([l1, r1]), max([l2, r2]), (vl + vr) / 2.0, (ampl + ampr) / 2.0, ...
                   phi, (dxl + dxr) / 2.0, (dyr + dyl) / 2.0, (displ + dispr) / 2.0];
            nsac = nsac + 1;
            labels = [labels; 1, 1];
        else
            labels = [labels; 1, 0];
            L_mono = [L_mono; L(i, :)];
        end
    end

    % Identify monocular microsaccades in right eye
    for i = 1:NR
        l1 = R(i, 1);
        l2 = R(i, 2);
        overlap = find((L(:, 2) >= l1) & (L(:, 1) <= l2));
        if isempty(overlap)
            labels = [labels; 0, 1];
            R_mono = [R_mono; R(i, :)];
        end
    end

    sac = sac;
    L_mono = L_mono;
    R_mono = R_mono;
    labels = labels;
end
