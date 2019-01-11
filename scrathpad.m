if all_combinations(1,this_comb_inside) == all_combinations(2,this_comb_inside)
    if this_clustsize > 1
        C=0;
        
        n=sum(neighbour_cells);
        cluster_membership{i}=zeros(n,1);
        
        visited=false(n,1);
        isnoise=false(n,1);
        
        %Calculate epsilon from the knee point of a
        %1000-bin histogram to the k-1th nearest neighbour
        %Remember that first index is self, so no need
        %to subtract one.
        [hist_y, hist_x] = hist(all_multi_real_distances(this_clustsize,:),1000);
        [epsilon, ~] = knee_pt(hist_y,hist_x,true);
        
        for this_cell=1:n
            if ~visited(this_cell)
                visited(this_cell)=true;
                
                disp(['Working on cell ' num2str(this_cell) ' out of ' num2str(n)]);
                Neighbors=setdiff(all_indexes(all_multi_real_distances(:,this_cell)<=epsilon,this_cell),this_cell)';
                if numel(Neighbors)<this_clustsize
                    % X(i,:) is NOISE
                    isnoise(this_cell)=true;
                else
                    C=C+1;
                    cluster_membership{i}(this_cell)=C;
                    k = 1;
                    while true
                        this_neighbour_cell = Neighbors(k);
                        
                        if ~visited(this_neighbour_cell)
                            visited(this_neighbour_cell)=true;
                            Neighbors2=setdiff(all_indexes(all_multi_real_distances(:,this_neighbour_cell)<=epsilon,this_neighbour_cell),this_neighbour_cell)';
                            if numel(Neighbors2)>=this_clustsize
                                Neighbors=[Neighbors Neighbors2];   %#ok
                            end
                        end
                        if cluster_membership{i}(this_neighbour_cell)==0
                            cluster_membership{i}(this_neighbour_cell)=C;
                        end
                        
                        k = k + 1;
                        if k > numel(Neighbors)
                            break;
                        end
                    end
                end
                
            end
        end
    end
end
%                     figure
%                     ax = gca();
%                     scatter(data{X_ind}(data{cell_ind}==1),data{Y_ind}(data{cell_ind}==1),1,cluster_membership{2}) %Ignore cell type 0
%                     hold(ax, 'on');
%                     imh = imshow(large_thumbnail_io);
%                     hold(ax, 'off');
%                     uistack(imh, 'bottom')
%
%                     test_cluster = cluster_membership{i};
%                     test_cluster(test_cluster~=0) = 1;
%
%                     figure
%                     ax = gca();
%                     scatter(data{X_ind}(data{cell_ind}==1),data{Y_ind}(data{cell_ind}==1),1,test_cluster) %Ignore cell type 0
%                     hold(ax, 'on');
%                     imh = imshow(large_thumbnail_io);
%                     hold(ax, 'off');
%                     uistack(imh, 'bottom')
%

x_data_here = data_trimmed{X_ind}(neighbour_cells);
y_data_here = data_trimmed{Y_ind}(neighbour_cells);

for this_cluster_number = 1:max(cluster_membership{i})
    this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
    this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
    this_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset);
    this_cluster_boundary{i}{this_cluster_number} = [this_x_subset(this_cluster_boundary_numbers),this_y_subset(this_cluster_boundary_numbers)];
end
figure
ax = gca();
hold on
for this_cluster_number = 1:length(this_cluster_boundary)
    if ~isempty(this_cluster_boundary{this_cluster_number})
        plot(this_cluster_boundary{this_cluster_number}(:,1),this_cluster_boundary{this_cluster_number}(:,2),'k-')
    end
end
hold(ax, 'on');
imh = imshow(large_thumbnail_io);
hold(ax, 'off');
uistack(imh, 'bottom')
title(['Minimum cluster size ' num2str(this_clustsize)]);


