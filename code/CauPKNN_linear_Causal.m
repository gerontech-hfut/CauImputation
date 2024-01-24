function imputed2 = CauPKNN_linear_Causal(missing,K_proportion)
complete = [];
no_complete = [];
no_complete_index = [];
complete_index = [];
[row, col] = size(missing);
Mi_dataset= zeros(col,col);
gauss_sum = zeros(col);
gauss_dataset=zeros(col,col);
logic_dataset=zeros(col,col);
total_K = 0;
total_count = 0; 

%extract complete data and incomplete data
for i = 1:row
    if isempty(find(missing(i,:)==9, 1))
        continue;
    else
        no_complete_index = [no_complete_index i];
    end
end

no_complete = missing(no_complete_index,:);
complete_index = setdiff(1:row,no_complete_index);
complete = missing(complete_index,:);

missing_rate = [];
for no_com = 1:length(no_complete_index)
    nan_index = find(no_complete(no_com,:)==9);
    nan_num = length(nan_index);
    nan_rate = nan_num/col;
    missing_rate = [missing_rate nan_rate];
end
[~,rate_rank] = sort(missing_rate);

% missing data imputation
for i = 1:length(rate_rank)
    if ~any(~(no_complete(rate_rank(i),:)==9))
        continue;
    end
    miss_row_index=rate_rank(i);
    miss_atr_index=find(no_complete(rate_rank(i),:)==9);
    miss_atr_num=length(miss_atr_index);
    %weight_array = [];
    column_index1=find(no_complete(rate_rank(i),:)~=9);%
    
    %op_column_index = ~column_index;
    for m = 1:miss_atr_num
        column_index = intersect(column_index1, find(Mi_dataset(miss_atr_index(m),:) ~= 0));
        temp = repmat(no_complete(rate_rank(i), column_index), length(complete_index), 1);
        distance = power((temp - complete(:, column_index)), 2);
        distance = distance * (Mi_dataset(miss_atr_index(m), column_index))';
        distance = sqrt(distance);
        %[~, distance_index] = sort(distance, 'descend');
        [~, distance_index] = sort(distance);
        make_K = sum(distance) * K_proportion;
        sum_distance = 0;
        K = 1;
        while(sum_distance < make_K)
            sum_distance = sum_distance + distance(distance_index(K));
            K = K + 1;
        end
        total_K = total_K + K;
        total_count = total_count + 1;
        adjusted_distance = distance(distance_index(1:K)) + 1e-10;

        %Modify the distance to make sure there is no zero value
        distance_p = 1 ./ adjusted_distance;
        distance_p = distance_p / sum(distance_p); 

        %Compute P( f_target = f_o);
        neighbor_values = complete(distance_index(1:K), miss_atr_index(m));
        unique_neighbor_values = unique(neighbor_values);
        value_counts = zeros(size(unique_neighbor_values));
        for uv = 1:length(unique_neighbor_values)
            value_counts(uv) = sum(complete(:, miss_atr_index(m)) == unique_neighbor_values(uv));
        end
        value_probs = value_counts / sum(value_counts);
        
        weighted_values = zeros(size(unique_neighbor_values));
        for ka = 1:K
            sample_value = neighbor_values(ka);
            value_index = find(unique_neighbor_values == sample_value);
            weighted_values(value_index) = weighted_values(value_index) + distance_p(ka) * length(unique_neighbor_values) * value_probs(value_index);
        end

        [~, max_index] = max(weighted_values);
        impute_value = unique_neighbor_values(max_index);
        

        %[~, max_index] = max(weighted_values);
        %impute_value = unique_neighbor_values(max_index);
        no_complete(rate_rank(i), miss_atr_index(m)) = impute_value;
    end
end
missing(no_complete_index,:) = no_complete;
imputed2 = missing;
% Calculate the average value of K
average_K = total_K / total_count;
fprintf('Average K: %f\n', average_K);