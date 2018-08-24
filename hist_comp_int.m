function hist_comp_int(base,neighbour)
% A function for plotting a histogram of observed and expected
% distributions between two cell types
% cell type options:
% key{1} = 'rubbish';
% key{2} = 'tumour';
% key{3} = 'lymphocyte';
% key{4} = 'stroma';
% key{5} = 'normal';

%Simulate output file header
all_combinations = combvec(0:4,0:4);
key{1} = 'rubbish';
key{2} = 'tumour';
key{3} = 'lymphocyte';
key{4} = 'stroma';
key{5} = 'normal';
header_string = [];
for this_comb = 1:size(all_combinations,2)
    header_string = [header_string ',Av_Mean_Distance_' key{all_combinations(1,this_comb)+1} '_to_' key{all_combinations(2,this_comb)+1} ',Av_Bootstrap_Distance_' key{all_combinations(1,this_comb)+1} '_to_' key{all_combinations(2,this_comb)+1} ',iqr_Mean_Distance_' key{all_combinations(1,this_comb)+1} '_to_' key{all_combinations(2,this_comb)+1} ',iqr_Bootstrap_Distance_' key{all_combinations(1,this_comb)+1} '_to_' key{all_combinations(2,this_comb)+1}];
end
full_string = ['Slide_ID,Cluster_Size,Num_Total,Num_Rubbish,Num_Tumour,Num_Lymphs,Num_Stroma,Num_Normal,Prop_Rubbish,Prop_Tumour,Prop_Lymphs,Prop_Stroma,Prop_Normal' header_string];
split_full_string = strsplit(full_string,',');

IndexC = strfind(split_full_string, ['Av_Mean_Distance_' base '_to_' neighbour]);
col_int = find(not(cellfun('isempty', IndexC)));
wei_data = csvread('clustering_data_nobootstrap.csv',1,0);

[n2,x2]=hist(wei_data(:,col_int),10000);
[n1,x1]=hist(wei_data(:,col_int+1)/100,10000);

figure
subplot(4,1,1)
h1=bar(x1,n1,'hist');
hold on
h2=bar(x2,n2,'hist');
set(h2,'facecolor','red')
set(h2,'edgecolor','none')
set(h1,'edgecolor','none')
alpha(0.3)
legend({'expected','observed'})
title(['Distance ' base ' to ' neighbour ])

subplot(4,1,2)
h1=bar(x1,n1,'hist');
hold on
h2=bar(x2,n2,'hist');
set(h2,'facecolor','red')
set(h2,'edgecolor','none')
set(h1,'edgecolor','none')
xlim([0 5000])
alpha(0.3)

subplot(4,1,3)
h1=bar(x1,n1,'hist');
hold on
h2=bar(x2,n2,'hist');
set(h2,'facecolor','red')
set(h2,'edgecolor','none')
set(h1,'edgecolor','none')
xlim([0 1000])
alpha(0.3)

subplot(4,1,4)
h1=bar(x1,n1,'hist');
hold on
h2=bar(x2,n2,'hist');
set(h2,'facecolor','red')
set(h2,'edgecolor','none')
set(h1,'edgecolor','none')
xlim([0 500])
alpha(0.3)

