function aggregate = aggregateData(directory)
    oldpath = addpath(['.',filesep,'MATLAB_scripts']);
    if nargin < 1
        directory = uigetdir('.','Select the directory containing the directories whose output should be aggregated');
    end

    disp(['Aggregating from ',directory,':'])

    [parentFolders, folders] = traverseDirs(directory, '*.tif_analysis*');

    aggregate = table();
    for i = 1 : length(parentFolders)
        try
            data = load([parentFolders{i}, filesep, folders{i},filesep,'output.mat']);
            data.parentDir = parentFolders{i};
        catch exception
            warning(['No output in ',parentFolders{i}, filesep, folders{i},'. Skipping.'])
            continue
        end
        % Find any non-singleton doubles and put them in a cell array.
        fn = fieldnames(data);
        for j = 1 : numel(fn)
            if (isnumeric(data.(fn{j})) || islogical(data.(fn{j}))) && numel(data.(fn{j})) > 1
                data.(fn{j}) = {data.(fn{j})};
            end
        end
        try 
            newTable = struct2table(data,'AsArray',true);
            newTable = [newTable(:,1), newTable(:,sort(newTable.Properties.VariableNames(2:end)))];
            newTable = assignDescriptionsAndUnitsToTable(newTable);
            newTable = removevars(newTable, toRemove());
            aggregate = [aggregate; newTable];
        catch exception
            warning(['Unable to append ',parentFolders{i}, filesep, folders{i},'. Skipping.'])
            continue;
        end
    end
    aggregate.Properties.RowNames = aggregate.filepath;
    disp([num2str(size(aggregate,1)),' files aggregated.'])

    % Sort the columns (variables) in alphabetical order.
    aggregate = [aggregate(:,1), aggregate(:,sort(aggregate.Properties.VariableNames(2:end)))];

    % disp(['Saving in ',directory,filesep,'aggregate.mat.',' and ',directory,filesep,'aggregate_summary.csv.'])
    % save([directory,filesep,'aggregate.mat'])
    writetable(aggregate,[directory,filesep,'aggregate.csv'],'Delimiter',',','QuoteStrings','all');
end