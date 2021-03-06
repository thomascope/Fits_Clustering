% A master script for analysing cell clustering in FITS data
% this_pool = cbupool(128);
% this_pool.ResourceTemplate = '-l nodes=^N^,mem=256GB,walltime=168:00:00';
% parpool(this_pool,this_pool.NumWorkers)

function Fits_Clustering_norubbish

cluster_size = [1, 2, 3, 4, 5, 10, 20, 30, 50]; % Define the size of the cluster of interest for distance in terms of number of cells
clusters_for_detail = [50]; %Define the sizes of the detailed outputs
only_detail = 1; %Don't process files unless they are in the list specified for detail;
if ~only_detail
    filenames = dir('./level2_catalogues/*.fits');
end

%detail_filenames = {'647365.fits','647366.fits','605012.fits','605019.fits','605181.fits','605182.fits','608225.fits','608226.fits','643632.fits','643619.fits','647364.fits','648121.fits'}; %Which files do we want detailed output for?

these_filenumbers = [
591564
591916
591535
591571
591541
591568
591877
591567
591843
591538
591554
590218
590263
590265
590274
590271
590268
592609
590822
590788
590785
590782
590779
590763
591874
612235
612091
612863
612862
612842
612836
612837
612828
612092
612093
612119
612120
612143
612144
612345
612328
612347
612348
616595
616575
616572
616573
616522
616843
616845
616848
616851
616854
616856
616859
616862
616871
616872
616892
616893
616976
616987
617482
602951
603006
603911
604105
604106
607166
626171
593978
593980
593981
593987
593988
593989
594000
594001
594002
594006
594007
594009
594017
594020
594022
594027
594033
594060
594061
594062
594065
594066
594068
594071
594072
594075
594080
594081
594082
594141
595806
595807
595808
597768
597769
597770
597777
597778
597779
597794
597795
597796
597800
597801
597802
599797
599798
599799
599803
599804
599806
599831
599832
599850
599851
599852
602915
602921
602935
602938
602942
602945
602952
602956
602958
602962
602966
602971
602976
602983
602994
603005
603007
603245
603257
603262
603263
603269
603271
603275
603919
603922
603923
603925
603928
603941
603944
603963
603966
603970
603976
603982
603987
604045
604047
604055
604057
604061
604062
604066
604071
604077
604078
604083
604089
604098
619460
619465
619483
619487
619503
619508
619515
619519
619527
619533
619550
619571
619572
619579
619842
619852
619859
619863
619868
619877
619884
619905
619922
619932
619942
619953
625322
625333
625338
625876
625882
625887
625891
625908
625911
625916
625923
625930
625931
625936
625946
625951
625958
626018
626102
626103
626160
626166
626172
650537
602927
597812
597813
597814
607165
593708
593942
593944
593949
593950
593952
593960
593961
593963
593971
593972
593974
594035
594036
594038
594044
594045
594046
594051
594052
594053
594088
594089
594091
594095
594096
594097
594105
594107
594110
594111
594113
594116
594117
594118
594137
594138
594140
597786
597788
597789
597806
597807
597808
597820
597821
597823
597830
597831
597832
597838
597839
598791
599813
599814
599817
599823
599824
599825
599839
599842
599843
599844
603055
603253
603279
603283
603288
603292
603298
603300
603356
603904
603905
603952
603956
603957
603958
604094
604107
604108
619470
619476
619489
619496
619538
619539
619540
619556
619558
619564
619565
619583
619678
619793
619799
619803
619809
619815
619820
619821
619832
619837
619847
619857
619872
619896
619897
619911
619916
619926
619937
619954
620961
625865
625871
625895
625900
626047
626061
650630
650631
673966
673975
673983
673996
674288
674276
674294
674305
674314
674327
674187
674188
673909
674197
674337
674215
674223
674231
674232
674118
674264
674130
674131
674344
674047
674143
674159
674160
684483
674056
];

detail_filenames = cell(size(these_filenumbers));
for i = 1:length(detail_filenames)
    detail_filenames{i} = [num2str(these_filenumbers(i)) '.fits'];
    if only_detail 
        filenames(i,1).name = detail_filenames{i};
    end
end


%Create output file
outfile = ['./clustering_data_multi_distance_third.csv'];
detail_dir = ['./detail_newlist/'];
if ~exist(detail_dir,'dir')
mkdir(detail_dir)
end
outfile_detail_stem = [detail_dir 'clustering_detail_'];
fileID = fopen(outfile,'w');

%Write output file header
all_combinations = combvec(1:4,1:4); % 0:4 includes rubbish, 1:4 excludes
key{1} = 'rubbish';
key{2} = 'tumour';
key{3} = 'lymphocyte';
key{4} = 'stroma';
key{5} = 'normal';
header_string = [];
for this_comb_outside = 1:size(all_combinations,2)
    header_string = [header_string ',Av_Mean_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1} ',Av_Bootstrap_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1} ',iqr_Mean_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1} ',iqr_Bootstrap_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1}];
end

fprintf(fileID,['Slide_ID,Cluster_Size,Num_Total,Num_Rubbish,Num_Tumour,Num_Lymphs,Num_Stroma,Num_Normal,Prop_Rubbish,Prop_Tumour,Prop_Lymphs,Prop_Stroma,Prop_Normal' header_string '\n']);
fclose(fileID);

parfor thisfile = 1:size(filenames,1) %could parallelise here by subject
%for thisfile = 1:size(filenames,1) %could parallelise here by subject
    %parfor thisfile = 1:5; %for testing
    if any(strcmp(filenames(thisfile).name,detail_filenames))
        detailed_output = 1;
    else
        detailed_output = 0;
    end
    try
        sprintf(['Working on file ' filenames(thisfile).name])
        
        this_comb_inside = 0
        
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
        
        % Now exclude anything marked as rubbish from the whole data
        data_trimmed = data;
        data_trimmed{cell_ind} = data_trimmed{cell_ind}(data{cell_ind}~=0);
        data_trimmed{X_ind} = data_trimmed{X_ind}(data{cell_ind}~=0);
        data_trimmed{Y_ind} = data_trimmed{Y_ind}(data{cell_ind}~=0);
        
        % Now compute the euclidian distances between each cell type and its
        % nearest neighbour of another cell type
        av_bootstrap_distance = zeros(length(cluster_size),size(all_combinations,2));
        iqr_bootstrap_distance = zeros(length(cluster_size),size(all_combinations,2));
        av_real_distance = zeros(length(cluster_size),size(all_combinations,2));
        iqr_real_distance = zeros(length(cluster_size),size(all_combinations,2));
        if detailed_output
            individual_real_distances = cell(length(clusters_for_detail),size(all_combinations,2));
        end
        
        for this_comb_inside = 1:size(all_combinations,2)
            sprintf(['Working on file ' filenames(thisfile).name ' combination ' num2str(this_comb_inside) ])
            base_cells = data_trimmed{cell_ind}==all_combinations(1,this_comb_inside);
            neighbour_cells = data_trimmed{cell_ind}==all_combinations(2,this_comb_inside);
            if sum(base_cells)==0||sum(neighbour_cells)==0
                av_bootstrap_distance(:,this_comb_inside) = NaN;
                av_real_distance(:,this_comb_inside) = NaN;
                iqr_bootstrap_distance(:,this_comb_inside) = NaN;
                iqr_real_distance(:,this_comb_inside) = NaN;
            else
                
                % First compute the real distances between the cell of interest and its
                % second nearest neighbour (to allow same cell type distance
                % computation)
                %av_real_distance(this_comb_inside) = mean(max(pdist2([data_trimmed{X_ind}(neighbour_cells) data_trimmed{Y_ind}(neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest',2)));
                %clear all_real_distances all_multi_real_distances
                [all_multi_real_distances] = pdist2([data_trimmed{X_ind}(neighbour_cells) data_trimmed{Y_ind}(neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest',max(cluster_size)+1);
                i = 0;
                j = 0;
                if size(all_multi_real_distances,1)<max(cluster_size)+1
                    all_multi_real_distances(size(all_multi_real_distances,1)+1:max(cluster_size)+1,:)=NaN;
                end
                for this_clustsize = cluster_size
                    i = i+1;
                    if any(all_multi_real_distances(1,:)==0) %Now exclude case where zero distance indicates it is the same point
                        all_mean_real_distances = mean(all_multi_real_distances(2:(this_clustsize+1),:),1);
                        all_real_distances = all_multi_real_distances(2:(this_clustsize+1),:);
                    else
                        all_mean_real_distances = mean(all_multi_real_distances(1:this_clustsize,:),1);
                        all_real_distances = all_multi_real_distances(1:(this_clustsize),:);
                    end
                    av_real_distance(i,this_comb_inside) = median(all_mean_real_distances);
                    iqr_real_distance(i,this_comb_inside) = iqr(all_mean_real_distances);
                    if detailed_output
                        if ismember(this_clustsize,clusters_for_detail)
                            j = j+1;
                            individual_real_distances{j,this_comb_inside}=all_real_distances;
                        end
                        
                        
                    end
                end
                
                
                % Then randomise the neighbours across all cell locations to create a
                % null distribution - just once here for computational
                % speed
                
                
                random_neighbour_cells = neighbour_cells(randperm(length(neighbour_cells)));
                [random_multi_distance] = pdist2([data_trimmed{X_ind}(random_neighbour_cells) data_trimmed{Y_ind}(random_neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest',max(cluster_size)+1);
                i = 0;
                if size(random_multi_distance,1)<max(cluster_size)+1
                    random_multi_distance(size(random_multi_distance,1)+1:max(cluster_size)+1,:)=NaN;
                end
                for this_clustsize = cluster_size
                    i = i+1;
                    if any(random_multi_distance(1,:)==0) %Now exclude case where zero distance indicates it is the same point
                        all_random_distances = mean(random_multi_distance(2:(this_clustsize+1),:),1);
                    else
                        all_random_distances = mean(random_multi_distance(1:this_clustsize,:),1);
                    end
                    av_bootstrap_distance(i,this_comb_inside) = median(all_random_distances);
                    iqr_bootstrap_distance(i,this_comb_inside) = iqr(all_random_distances);
                end
                
            end
        end
        i = 0;
        outfile = ['./clustering_data_multi_distance_third.csv'];
        fileID = fopen(outfile,'a');
        for this_clustsize = cluster_size
            i = i+1;
            data_string = [];
            for this_comb_inside = 1:size(all_combinations,2)
                data_string = [data_string ',' num2str(av_real_distance(i,this_comb_inside)) ',' num2str(av_bootstrap_distance(i,this_comb_inside)) ',' num2str(iqr_real_distance(i,this_comb_inside)) ',' num2str(iqr_bootstrap_distance(i,this_comb_inside))];
            end
            
            fprintf(fileID,[filenames(thisfile).name(1:end-5) ',' num2str(cluster_size(i)) ',' num2str(num_total) ',' num2str(num_rubbish) ',' num2str(num_tum_cells) ',' num2str(num_ly_cells) ',' num2str(num_str_cells) ',' num2str(num_norm_cells) ',' num2str(prop_rubbish) ',' num2str(prop_tum_cells) ',' num2str(prop_ly_cells) ',' num2str(prop_str_cells) ',' num2str(prop_norm_cells) data_string '\n']);
        end
        fclose(fileID);
        
        if detailed_output
            for this_comb_inside = 1:size(all_combinations,2)
                source_slide = strsplit(filenames(thisfile).name,'.fits');
                %header_string = [header_string ',Av_Mean_Distance_' key{all_combinations(1,this_comb_inside)+1} '_to_' key{all_combinations(2,this_comb_inside)+1} ',Av_Bootstrap_Distance_' key{all_combinations(1,this_comb_inside)+1} '_to_' key{all_combinations(2,this_comb_inside)+1} ',iqr_Mean_Distance_' key{all_combinations(1,this_comb_inside)+1} '_to_' key{all_combinations(2,this_comb_inside)+1} ',iqr_Bootstrap_Distance_' key{all_combinations(1,this_comb_inside)+1} '_to_' key{all_combinations(2,this_comb_inside)+1}];
                j = 0;
                for this_clustsize = clusters_for_detail
                    j = j+1;
                    this_fname = [outfile_detail_stem source_slide{1} '_' key{all_combinations(1,this_comb_inside)+1} '_to_' key{all_combinations(2,this_comb_inside)+1} '_nearest' num2str(clusters_for_detail(j)) '.csv'];

                    csvwrite(this_fname,individual_real_distances{j,this_comb_inside}');
                    
                    sprintf(['Working on file ' filenames(thisfile).name ' nearest ' num2str(this_clustsize(j)) ])
                    [all_multi_real_distances, all_multi_real_indexes] = pdist2([data_trimmed{X_ind} data_trimmed{Y_ind}],[data_trimmed{X_ind} data_trimmed{Y_ind}],'euclidean','Smallest',max(cluster_size)+1);
                    if any(all_multi_real_distances(1,:)==0) %Now exclude case where zero distance indicates it is the same point
                        all_real_indexes = all_multi_real_indexes(2:(this_clustsize+1),:);
                        all_real_distances = all_multi_real_distances(2:(this_clustsize+1),:);
                    else
                        all_real_indexes = all_multi_real_indexes(1:this_clustsize,:);
                        all_real_distances = all_multi_real_distances(1:this_clustsize,:);
                    end
                    
                    for this_cell_type = 1:4
                        this_fname = [outfile_detail_stem source_slide{1} '_' key{this_cell_type+1} '_cellids_nearest' num2str(clusters_for_detail(j)) '.csv'];
                        these_real_indexes = all_real_indexes(:,data_trimmed{cell_ind}==this_cell_type);
                        cell_ids_out = nan(size(these_real_indexes));
                        for this_cell = 1:this_clustsize
                            cell_ids_out(this_cell,:)=data_trimmed{cell_ind}(these_real_indexes(this_cell,:));
                        end
                        csvwrite(this_fname,cell_ids_out');
                        
                        this_fname = [outfile_detail_stem source_slide{1} '_' key{this_cell_type+1} '_distances_nearest' num2str(clusters_for_detail(j)) '.csv'];
                        these_real_distances = all_real_distances(:,data_trimmed{cell_ind}==this_cell_type);
                        csvwrite(this_fname,these_real_distances');
                    end
                                       
                    
                end
            end
                
            
        end
        sprintf(['Finished file ' filenames(thisfile).name])
        
    catch
        outfile = ['./clustering_data_multi_distance_third.csv'];
        fileID = fopen(outfile,'a');
        if this_comb_inside > 0
            fprintf(fileID,[filenames(thisfile).name(1:end-5) ',failed at ' num2str(this_comb_inside) '\n']);
        else
            fprintf(fileID,[filenames(thisfile).name(1:end-5) ' does not exist \n']);
        end
        fclose(fileID);
        sprintf(['Failed file ' filenames(thisfile).name ' moving on'])
    end
end

%end
