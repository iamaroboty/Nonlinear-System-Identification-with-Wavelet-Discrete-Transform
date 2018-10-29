%%%%%%%%%%%%%%  StructVsCell.m %%%%%%%%%%%%%%%

clear all

M = 100; % number of repetitions
N = 2^10; % size of cell array and struct


for m = 1:M
    % Fill up a template cell array with
    % lists of randomly sized matrices with
    % random elements.
    template{N} = 0;
    for n = 1:N
        r1 = round(24*rand());
        r2 = round(24*rand());
        r3 = rand(round(r2*rand),round(r1*rand()));
        template{N} = r3;
    end

    % Make a cell array equivalent
    % to the template.
    cell_array = template;

    % Create a struct with the
    % same data.
    structure = struct('data',0);
    for n = 1:N
        structure(n).data = template{n};
    end

    % Time cell array
    tic;
    for n = 1:N
        data = cell_array{n};
        cell_array{n} = data';
    end
    cell_time(m) = toc;

    % Time struct
    tic;
    for n = 1:N
        data = structure(n).data;
        structure(n).data = data';
    end
    struct_time(m) = toc;
end

str = sprintf('cell array: %0.4f',mean(cell_time));
disp(str);
str = sprintf('struct: %0.4f',mean(struct_time));
disp(str);
str = sprintf('struct_time / cell_time: %0.4f',mean(struct_time)/mean(cell_time));
disp(str);

% Check memory use
whos
