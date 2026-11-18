% Patient  Demographic covariates shared by every layer.
%
%   This file is reproduced from the host simulator so that the
%   supplementary package is self-contained when distributed on its own.
%   setupPaths adds this folder only when Patient is not already on the
%   MATLAB path, so the host definition takes precedence inside the
%   parent repository.

classdef Patient
    properties (SetAccess = immutable)
        id
        age     % years
        height  % cm
        weight  % kg
        sex     % M - male, F - female
        bmi
        lbm
    end

    methods
        function obj = Patient(id, age, height, weight, sex)
            obj.id = id;
            obj.age = age;
            obj.height = height;
            obj.weight = weight;
            obj.sex = sex;
            obj.bmi = weight / ((height/100)^2);

            if iscell(sex)
                sexChar = sex{1};
            else
                sexChar = char(sex);
            end

            switch sexChar 
                case 'M'
                    obj.lbm = 1.1 * weight - 128 * (weight/height)^2; % James Formula for Men
                case 'F'
                    obj.lbm = 1.07 * weight - 148 * (weight/height)^2; % James Formula for Women
                otherwise
                    disp('Warning: Undefined gender, lbm computed with male formula');
                    obj.lbm = 1.1 * weight - 128 * (weight/height)^2; % James Formula for Men
            end
        end
    end

    methods (Static)
        function obj = fromVitalDbStruct(s)
            if numel(s) > 1
                obj(numel(s),1) = Patient(0,0,0,0,'M'); % Preallocate
                for i = 1:numel(s)
                    obj(i) = Patient.fromVitalDbStruct(s(i));
                end
                return
            end

            requiredFields = {'caseid','age','height','weight','sex'};
            for f = requiredFields
                if ~isfield(s, f{1})
                    error('Field "%s" missing in input struct.', f{1});
                end
            end

            obj = Patient(...
                s.caseid, ...
                s.age, ...
                s.height, ...
                s.weight, ...
                s.sex ...
            );
        end
    end

end