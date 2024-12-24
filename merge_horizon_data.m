function [all_data, subj_mapping] = merge_horizon_data(root, ids, groupdata, load_type)        
    all_data = cell(1, numel(ids)); 
    
    good_index = [];
    subj_mapping = cell(numel(ids)*numel(load_type), 2); 
    folder = './data';
    n=0;
    for i = 1:numel(ids)
        for j = 1:numel(load_type)
        n = n+1;
        load=load_type{j};
        id = ids{i};  

        all_data{n} = compile_data(folder, id, groupdata, load);  
        
        if (ismember(size(all_data{n}, 1), 80) && (ismember(sum(all_data{n}.gameLength), 600)))
            good_index = [good_index n];
        end
        
        all_data{n}.subjectID = repmat(n, size(all_data{n}, 1), 1);
        
        subj_mapping{n, 1} = {[id '-' num2str(j)]};
        subj_mapping{n, 2} = n;
        end
    end
    
    all_data = all_data(good_index);
    subj_mapping = subj_mapping(good_index, :);
    
    all_data = vertcat(all_data{:});    
end