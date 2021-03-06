% A master script for analysing cell clustering in FITS data
%cbupool(128) % Open parallel worker pool

filenames = dir('./level2_catalogues/*.fits');
cluster_size = 5; % Define the size of the cluster of interest for distance in terms of number of cells

%Create output file
outfile = ['./clustering_data.csv'];
fileID = fopen(outfile,'w');

%Write output file header
all_combinations = combvec(0:4,0:4);
key{1} = 'rubbish';
key{2} = 'tumour';
key{3} = 'lymphocyte';
key{4} = 'stroma';
key{5} = 'normal';
header_string = [];
for this_comb = 1:size(all_combinations,2)
    header_string = [header_string ',Av_Mean_Distance_' key{all_combinations(1,this_comb)+1} '_to_' key{all_combinations(2,this_comb)+1} ',Percentile_Mean_Distance_' key{all_combinations(1,this_comb)+1} '_to_' key{all_combinations(2,this_comb)+1} ',Zscore_Mean_Distance_' key{all_combinations(1,this_comb)+1} '_to_' key{all_combinations(2,this_comb)+1}];
end

fprintf(fileID,['Slide_ID,Cluster_Size,Num_Total,Num_Rubbish,Num_Tumour,Num_Lymphs,Num_Stroma,Num_Normal,Prop_Rubbish,Prop_Tumour,Prop_Lymphs,Prop_Stroma,Prop_Normal' header_string '\n']);
fclose(fileID);


for thisfile = 1:size(filenames,1) %could parallelise here by subject
    %parfor thisfile = 1:5; %for testing
    try
        sprintf(['Working on file ' filenames(thisfile).name])
        
        data = fitsread(['./level2_catalogues/' filenames(thisfile).name],'binarytable');
        info = fitsinfo(['./level2_catalogues/' filenames(thisfile).name]);
        
        % Assume that X_global is always the third field, Y_global the fourth
        % field, and cell type the final field.
        
        X_ind = 3;
        Y_ind = 4;
        cell_ind = size(data,2);
        
        %assert(max(data{cell_ind})==4&&min(data{cell_ind})==0,['Cell types must span 0-4 and they do not in file ' filenames(thisfile).name])
        % 0 is rubbish
        % 1 is tumour
        % 2 is lymphocyte
        % 3 is stroma
        % 4 is normal
        
        % %Optionally visualise the slide
        % figure
        % scatter(data{X_ind}(data{cell_ind}~=0),data{Y_ind}(data{cell_ind}~=0),1,data{cell_ind}(data{cell_ind}~=0)) %Ignore cell type 0
        
        % First compute the number of each cell type
        num_total = size(data{cell_ind},1);
        num_rubbish = sum(data{cell_ind}==0);
        num_tum_cells = sum(data{cell_ind}==1);
        num_ly_cells = sum(data{cell_ind}==2);
        num_str_cells = sum(data{cell_ind}==3);
        num_norm_cells = sum(data{cell_ind}==4);
        
        % Then compute the proportion of each cell type
        prop_rubbish = num_rubbish/num_total;
        prop_tum_cells = num_tum_cells/num_total;
        prop_ly_cells = num_ly_cells/num_total;
        prop_str_cells = num_str_cells/num_total;
        prop_norm_cells = num_norm_cells/num_total;
        
        % Now compute the euclidian distances between each cell type and its
        % nearest neighbour of another cell type
        distance_statsig = zeros(1,size(all_combinations,2));
        av_real_distance = zeros(1,size(all_combinations,2));
        for this_comb = 1:size(all_combinations,2)
            sprintf(['Working on file ' filenames(thisfile).name ' combination ' num2str(this_comb) ])
            base_cells = data{cell_ind}==all_combinations(1,this_comb);
            neighbour_cells = data{cell_ind}==all_combinations(2,this_comb);
            if sum(base_cells)==0||sum(neighbour_cells)==0
                distance_statsig(this_comb) = NaN;
                av_real_distance(this_comb) = NaN;
                
            else
                
                % First compute the real distances between the cell of interest and its
                % second nearest neighbour (to allow same cell type distance
                % computation)
                %av_real_distance(this_comb) = mean(max(pdist2([data{X_ind}(neighbour_cells) data{Y_ind}(neighbour_cells)],[data{X_ind}(base_cells) data{Y_ind}(base_cells)],'euclidean','Smallest',2)));
                clear all_real_distances all_multi_real_distances
                all_multi_real_distances(:,:) = pdist2([data{X_ind}(neighbour_cells) data{Y_ind}(neighbour_cells)],[data{X_ind}(base_cells) data{Y_ind}(base_cells)],'euclidean','Smallest',cluster_size+1);
                if any(all_multi_real_distances(1,:)==0) %Now exclude case where zero distance indicates it is the same point
                    all_real_distances = mean(all_multi_real_distances(2:end,:),1);
                else
                    all_real_distances = mean(all_multi_real_distances(1:end-1,:),1);
                end
                av_real_distance(this_comb) = median(all_real_distances);
                
                % Then randomise the neighbours across all cell locations to create a
                % null distribution
                num_bootstraps = 100; % Increase to 1000 for real analysis
                random_distance = zeros(1,num_bootstraps);
                parfor random_perm = 1:num_bootstraps % Or parallelise here by bootstrap
                    random_neighbour_cells = neighbour_cells(randperm(length(neighbour_cells)));
                    random_multi_distance = pdist2([data{X_ind}(random_neighbour_cells) data{Y_ind}(random_neighbour_cells)],[data{X_ind}(base_cells) data{Y_ind}(base_cells)],'euclidean','Smallest',cluster_size+1);
                    if any(random_multi_distance(1,:)==0) %Now exclude case where zero distance indicates it is the same point
                        all_random_distances = mean(random_multi_distance(2:end,:),1);
                    else
                        all_random_distances = mean(random_multi_distance(1:end-1,:),1);
                    end
                    random_distance(random_perm) = median(all_random_distances);
                end
                
                %Now compute the relative rank of the observed data in the randomised
                %data
                this_dataset=[av_real_distance(this_comb),random_distance];
                
                distance_statsig(this_comb)=sum(this_dataset(:)<av_real_distance(this_comb))/numel(this_dataset);
                
                distance_zscore(this_comb)=(av_real_distance(this_comb)-mean(random_distance))/std(random_distance);
                
            end
        end
        
        data_string = [];
        for this_comb = 1:size(all_combinations,2)
            data_string = [data_string ',' num2str(av_real_distance(this_comb)) ',' num2str(distance_statsig(this_comb)*100) ',' num2str(distance_zscore(this_comb))];
        end
        
        outfile = ['./clustering_data.csv'];
        fileID = fopen(outfile,'a');
        fprintf(fileID,[filenames(thisfile).name(1:end-5) ',' num2str(cluster_size) ',' num2str(num_total) ',' num2str(num_rubbish) ',' num2str(num_tum_cells) ',' num2str(num_ly_cells) ',' num2str(num_str_cells) ',' num2str(num_norm_cells) ',' num2str(prop_rubbish) ',' num2str(prop_tum_cells) ',' num2str(prop_ly_cells) ',' num2str(prop_str_cells) ',' num2str(prop_norm_cells) data_string '\n']);
        fclose(fileID);
        sprintf(['Finished file ' filenames(thisfile).name])
    catch
        outfile = ['./clustering_data.csv'];
        fileID = fopen(outfile,'a');
        fprintf(fileID,[filenames(thisfile).name(1:end-5) ',failed at ' num2str(this_comb) '\n']);
        fclose(fileID);
        sprintf(['Failed file ' filenames(thisfile).name ' moving on'])
    end
end

%end
