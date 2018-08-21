% A master script for analysing cell clustering in FITS data

filenames = dir('./level2_catalogues/*.fits');

%Create output file
outfile = ['./clustering_data.csv'];
fileID = fopen(outfile,'w');

%Write output file header
fprintf(fileID,['Slide_ID,Num_Total,Num_Rubbish,Num_Tumour,Num_Lymphs,Num_Stroma,Num_Normal,Prop_Rubbish,Prop_Tumour,Prop_Lymphs,Prop_Stroma,Prop_Normal\n']);
fclose(fileID);
%cbupool(64) % Open parallel worker pool

%parfor thisfile = 1:size(filenames,1) %for parallel
parfor thisfile = 1:5; %for testing

data = fitsread(['./level2_catalogues/' filenames(thisfile).name],'binarytable');
info = fitsinfo(['./level2_catalogues/' filenames(thisfile).name]);

% Assume that X_global is always the third field, Y_global the fourth
% field, and cell type the final field.

X_ind = 3;
Y_ind = 4;
cell_ind = size(data,2);

assert(max(data{cell_ind})==4&&min(data{cell_ind})==0,['Cell types must span 0-4 and they do not in file ' filenames(thisfile).name])
% 0 is rubbish
% 1 is tumour
% 2 is lymphocyte
% 3 is stroma
% 4 is normal

% %Optionally visualise the slide
% figure
% scatter(data{X_ind}(data{cell_ind}~=0),data{Y_ind}(data{cell_ind}~=0),1,data{cell_ind}(data{cell_ind}~=0)) %Ignore cell type 0

% First compute the number of each cell type
num_total = size(data{cell_ind},1)
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
all_combinations = combvec(0:4,0:4);
distance_statsig = zeros(1,size(all_combinations,2));
for this_comb = 1:size(all_combinations,2)
    
    base_cells = data{cell_ind}==all_combinations(1,this_comb);
    neighbour_cells = data{cell_ind}==all_combinations(2,this_comb);
    
    % First compute the real distances between the cell of interest and its
    % nearest neighbour
    real_distances = pdist2([data{X_ind}(neighbour_cells) data{Y_ind}(neighbour_cells)],[data{X_ind}(base_cells) data{Y_ind}(base_cells)],'euclidean','Smallest',1);
    
    % Then randomise the neighbours across all cell locations to create a
    % null distribution
    num_bootstraps = 100; % Increase to 1000 for real analysis
    random_distances = zeros(num_bootstraps,size(real_distances,2));
    for random_perm = 1:num_bootstraps
        random_distances(random_perm,:) = pdist2([data{X_ind}(neighbour_cells(randperm(length(neighbour_cells)))) data{Y_ind}(neighbour_cells(randperm(length(neighbour_cells))))],[data{X_ind}(base_cells) data{Y_ind}(base_cells)],'euclidean','Smallest',1);
    end
    
    %Now compute the relative rank of the observed data in the randomised
    %data
    this_dataset=[mean(real_distances);mean(random_distances,2)];
    
    distance_statsig(this_comb)=sum(this_dataset(:)<mean(real_distances))/numel(this_dataset);
    
    
end




outfile = ['./clustering_data.csv'];
fileID = fopen(outfile,'a');
fprintf(fileID,[filenames(thisfile).name(1:end-5) ',' num2str(num_total) ',' num2str(num_rubbish) ',' num2str(num_tum_cells) ',' num2str(num_ly_cells) ',' num2str(num_str_cells) ',' num2str(num_norm_cells) ',' num2str(prop_rubbish) ',' num2str(prop_tum_cells) ',' num2str(prop_ly_cells) ',' num2str(prop_str_cells) ',' num2str(prop_norm_cells) '\n']);
fclose(fileID);
end

%end
