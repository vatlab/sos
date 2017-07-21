% This file is part of Script of Scripts (sos), a workflow system
% for the execution of commands and scripts in different languages.
% Please visit https://github.com/vatlab/SOS for more information.
%
% Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

function [repr] = sos_py_repr (obj)
    if isnumeric(obj)

        if isinf(obj)
            if obj>0
                repr='np.inf';
            else
                repr='-np.inf';
            end

        elseif isreal(obj)==0
            rl=num2str(real(obj));
            im=num2str(imag(obj));
            repr=strcat('complex(',rl,',',im,')');

        elseif isnan(obj)
            repr='None';

        elseif isvector(obj)
            if length(obj)==1
                repr=num2str(obj);

            else
                nu=sprintf('%.0f,' ,obj);
                nu=nu(1:end-1);
                repr=strcat('np.array([',nu,'])');
            end

        elseif ismatrix(obj)
            save(fullfile(tempdir,'mat2py.mat'),'obj');
            %dirpath=tempdir;
            %repr="mat_contents = sio.loadmat('/Users/manchongleong/Desktop/mat2py.mat')"
            repr=char(strcat("sio.loadmat('",tempdir,"mat2py.mat')","['","obj","']"));

        else
            repr=num2str(obj);
        end

    elseif isstruct(obj)
        save(fullfile(tempdir,'stru2py.mat'),'obj');
        dirpath=tempdir;
        %repr="stru2py = sio.loadmat('/Users/manchongleong/Desktop/stru2py.mat')";
        repr=strcat("sio.loadmat('",tempdir,"stru2py.mat')","['","obj","']");

    elseif iscell(obj)
        save(fullfile(tempdir,'cell2py.mat'),'obj');
        dirpath=tempdir;
        %repr="stru2py = sio.loadmat('/Users/manchongleong/Desktop/stru2py.mat')";
        repr=strcat("sio.loadmat('",tempdir,"cell2py.mat')","['","obj","']");

    elseif islogical(obj)
        if obj
            repr = 'True';
        else
            repr = 'False';
        end

    else
        % unrecognized/unsupported datatype is transferred from
        % matlab to Python as None
        repr = 'None';
    end
end

