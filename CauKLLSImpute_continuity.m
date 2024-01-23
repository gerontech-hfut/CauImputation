function imputed2 = CauKLLSImpute_continuity(missing,sigma,K_proportion)
complete = [];
no_complete = [];
no_complete_index = [];
complete_index = [];
[row, col] = size(missing);
Mi_dataset= zeros(col,col);
gauss_sum = zeros(col);
gauss_dataset=zeros(col,col);
%logic_dataset=zeros(col,col);   

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


% Calculate the correlation matrix
for k=1:col
    sum_i=0;
    for j=1:col
        if k~=j
            sum_i=sum_i + abs(corr(missing(:,j),missing(:,k),'type', 'Kendall'));
        end
    end
    for j=1:col
        if k~=j
            Mi_dataset(k,j)=abs(corr(missing(:,j),missing(:,k),'type', 'Kendall')/sum_i);
        else
            Mi_dataset(k, j) = max(Mi_dataset(j,:)); 
        end
    end
end
%Best RMSE: 0.271504, Best sigmaor: 0.040000,Best rateor: 1.000000 Best K_proportionor: 0.200000，Best sigma_2: 0.250000

for k=1:col
    gauss_sum(k)=0;
    for j=1:col
        gauss_dataset(k,j)=RBF_function(1,Mi_dataset(k,j),sigma);
        gauss_sum(k)=gauss_sum(k)+gauss_dataset(k,j);
    end
end

gauss_sum = round((gauss_sum/max(gauss_sum))*col+col);
gauss_sum(gauss_sum>col)=col;
Mi_dataset= zeros(col,col);

for k=1:col    
    real_index=MAX_K_NUMBER(gauss_dataset(k,:),gauss_sum(k));
    Mi_dataset(k,real_index)=gauss_dataset(k,real_index);
end
    
missing_rate = [];
for no_com = 1:length(no_complete_index)
    nan_index = find(no_complete(no_com,:)==9);
    nan_num = length(nan_index);
    nan_rate = nan_num/col;
end
[~,rate_rank] = sort(missing_rate);
SIGMA_DE2 = 2;

% missing data imputation
for i = 1:length(rate_rank)
    if ~any(~(no_complete(rate_rank(i),:)==9))
        continue;
    end
    miss_row_index=rate_rank(i);
    miss_atr_index=find(no_complete(rate_rank(i),:)==9);
    miss_atr_num=length(miss_atr_index);
    %weight_array = [];
    column_index1=find(no_complete(rate_rank(i),:)~=9);
    
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
        while sum_distance < make_K && K <= length(distance)
            sum_distance = sum_distance + distance(distance_index(K));
            K = K + 1;
        end

        neighbor_values = complete(distance_index(1:K), miss_atr_index(m));
        
        % Modify the distance to make sure there is no zero value
        adjusted_distance = distance(distance_index(1:K)) + 1e-10;
        distance_p = 1 ./ adjusted_distance;
        distance_p = distance_p / sum(distance_p); 
    
        % Discrete Continuous Values
        %num_bins = K;
        num_bins = ceil(sqrt(row));
        bin_edges = linspace(min(complete(:, miss_atr_index(m))), max(complete(:, miss_atr_index(m))), num_bins + 1);
        binned_neighbor_values = discretize(neighbor_values, bin_edges);
    
        % Compute P( f_target = f_o);
        value_counts = histcounts(complete(:, miss_atr_index(m)), bin_edges);
        value_probs = value_counts / sum(value_counts);
    
        % constructing the weighting matrix CW
        weight_matrix = zeros(K, K);
        for ka = 1:K
            bin_index = binned_neighbor_values(ka);
            if ~isnan(bin_index)
                weight_matrix(ka, ka) = distance_p(ka) * value_probs(bin_index);
            end
        end
        for ka = 1:K
            weight_matrix(ka, ka) = RBF_function(weight_matrix(ka, ka), 0, SIGMA_DE2);
        end

        pearson_A = complete(distance_index(1:K), :);
        x = pinv((pearson_A(:, column_index)' * weight_matrix)) * (no_complete(rate_rank(i), column_index)');
        predict_value = x' * weight_matrix * pearson_A(:, miss_atr_index(m));
        no_complete(rate_rank(i), miss_atr_index(m)) = predict_value;
    end
end
missing(no_complete_index,:) = no_complete;
imputed2 = missing;