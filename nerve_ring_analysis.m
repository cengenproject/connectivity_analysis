clear all
clc
close all

addpath(genpath('./'))

%% USER OPTIONS
gene_family = 'cams'; %the gene family to investigate, options: 'cams', 'innexins', 'homeobox', 'gpcrs', 'nhrs'
neuron_1 = 'RIP'; %first neuron to investigate
neuron_2 = []; %second neuron to inverstigate (if blank [], will perform neuron 1 + synaptic partners vs. neuron 1 + adjacent partners test
pre_synaptic = 1; %If 1, will interrogate pre-synaptic partners vs. neuron, else will look at post-synaptic partners
include_self_gene_expr = 1; %If 1, will count the neuron 1/2's own gene expression level in the synaptic partnership gene count
unique_partners = 0; %If 1, for the neuron vs. neuron comparison, will only compara unique synaptic partners of neurons -- this doesn't apply to neuron synaptic partners vs. neuron adjacent partners since they are mutually exclusive anyway
number_of_heterophilic_genes= 20; %number of heterophilic interactions to plot for each group
%% INPUT FILES

column_gene_expr_file='./gene_expr_data/070120_Conn_matching_no_over_expression_averaged_neuron.csv';
row_gene_expr_file='./gene_expr_data/070120_Conn_matching_no_over_expression_averaged_neuron.csv';

chemical_column_neurons_file='./connectome_data/102820_nerve_ring_connectome_column_neurons.csv';
chemical_row_neurons_file='./connectome_data/102820_nerve_ring_connectome_row_neurons.csv';
chemical_connectivity_matrix_file='./connectome_data/102820_nerve_ring_connectome.csv';

adjacency_column_neurons_file='./connectome_data/102820_nerve_ring_adjacency_column_neurons.csv';
adjacency_row_neurons_file='./connectome_data/102820_nerve_ring_adjacency_row_neurons.csv';
adjacency_matrix_file='./connectome_data/102820_nerve_ring_adjacency.csv';

neuron_class_size_file='./gene_expr_data/Number_of_neurons_per_single_cell_class.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('importing csv files...')
column_expr_table=readtable(column_gene_expr_file);
column_expr_neurons=column_expr_table.Properties.VariableNames(2:end);
column_expr_genes=table2array(column_expr_table(:,1));
column_expr_matrix=table2array(column_expr_table(:,2:end));

row_expr_table=readtable(row_gene_expr_file);
row_expr_neurons=row_expr_table.Properties.VariableNames(2:end);
row_expr_genes=table2array(row_expr_table(:,1));
row_expr_matrix=table2array(row_expr_table(:,2:end));


adjacency_column_neurons=readtable(adjacency_column_neurons_file,'readvariablenames',false);
adjacency_conn_matrix_column_neurons=table2array(adjacency_column_neurons);
adjacency_row_neurons=readtable(adjacency_row_neurons_file,'readvariablenames',false);
adjacency_conn_matrix_row_neurons=table2array(adjacency_row_neurons);
adjacency_matrix_table=csvread(adjacency_matrix_file);

chemical_column_neurons=readtable(chemical_column_neurons_file,'readvariablenames',false);
chemical_conn_matrix_column_neurons=table2array(chemical_column_neurons);
chemical_row_neurons=readtable(chemical_row_neurons_file,'readvariablenames',false);
chemical_conn_matrix_row_neurons=table2array(chemical_row_neurons);
chemical_conn_matrix_table=csvread(chemical_connectivity_matrix_file);

neuron_class_size_table = readtable(neuron_class_size_file);
neuron_class_size_neurons=table2array(neuron_class_size_table(:,1));
neuron_class_size=table2array(neuron_class_size_table(:,2));


disp('importing csv files (done)')

disp('making sure connectivity matrix and gene expr. neurons match up...')


[a,b]=ismember(adjacency_conn_matrix_row_neurons,adjacency_conn_matrix_column_neurons);

adjacency_conn_matrix_column_neurons=adjacency_conn_matrix_column_neurons(b(a==1));
adjacency_matrix_table=adjacency_matrix_table(:,b(a==1));

[a,b]=ismember(adjacency_conn_matrix_row_neurons,row_expr_neurons);
[c,d]=ismember(adjacency_conn_matrix_row_neurons,chemical_conn_matrix_row_neurons);
['gene expr. row neurons' ' ' 'conn. matrix row neurons']
[row_expr_neurons(b(a==1))' chemical_conn_matrix_row_neurons(d(c==1)) adjacency_conn_matrix_row_neurons]


row_expr_matrix=row_expr_matrix(:,b(a==1));
row_expr_neurons=row_expr_neurons(b(a==1));

chemical_conn_matrix_table=chemical_conn_matrix_table(d(c==1),:);
chemical_conn_matrix_row_neurons=chemical_conn_matrix_row_neurons(d(c==1));


[a,b]=ismember(adjacency_conn_matrix_column_neurons,column_expr_neurons);
[c,d]=ismember(adjacency_conn_matrix_column_neurons,chemical_conn_matrix_column_neurons);
['gene expr. column neurons' ' ' 'conn. matrix column neurons']
[column_expr_neurons(b(a==1))' chemical_conn_matrix_column_neurons(d(c==1)) adjacency_conn_matrix_column_neurons]


column_expr_matrix=column_expr_matrix(:,b(a==1));
column_expr_neurons=column_expr_neurons(b(a==1));


chemical_conn_matrix_table=chemical_conn_matrix_table(:,d(c==1));
chemical_conn_matrix_column_neurons=chemical_conn_matrix_column_neurons(d(c==1));



[a,b]=ismember(adjacency_conn_matrix_row_neurons,neuron_class_size_neurons);

neuron_class_size_neurons=neuron_class_size_neurons(b(a==1));
neuron_class_size=neuron_class_size(b(a==1));

[neuron_class_size_neurons adjacency_conn_matrix_column_neurons]

C_adjacency=adjacency_matrix_table;
C_chemical=chemical_conn_matrix_table;
C_adjacency = C_adjacency./(neuron_class_size*neuron_class_size');
C_chemical = C_chemical./(neuron_class_size*neuron_class_size');

G1=row_expr_matrix';
G2=column_expr_matrix';

[a,b]=ismember(adjacency_conn_matrix_row_neurons,neuron_class_size_neurons);
['class size neurons' ' ' 'conn. matrix column neurons']
[neuron_class_size_neurons(b(a==1)) adjacency_conn_matrix_row_neurons]
neuron_class_size_neurons=neuron_class_size_neurons(b(a==1));
neuron_class_size=neuron_class_size(b(a==1));





disp(['Connectivity matrix: ' num2str(size(C_adjacency,1)) 'x' num2str(size(C_adjacency,2))])
disp(['Row-gene matrix: ' num2str(size(G1,1)) 'x' num2str(size(G1,2))])
disp(['Column-gene matrix: ' num2str(size(G2,1)) 'x' num2str(size(G2,2))])


load innexins.mat
load nhrs.mat
load cehs.mat
load homeobox.mat
load gpcrs.mat
load cams.mat

if strcmpi(gene_family,'cams')
    gene_set=cams;
elseif strcmpi(gene_family,'innexins')
    gene_set=innexins;
elseif strcmpi(gene_family,'homeobox')
    gene_set=homeobox;
elseif strcmpi(gene_family,'gpcrs')
    gene_set=gpcrs;
elseif strcmpi(gene_family,'nhrs')
    gene_set=nhrs;
end

[a,b]=ismember(gene_set,row_expr_genes);
G1=G1(:,b(a==1));
row_expr_genes=row_expr_genes(b(a==1));
[a,b]=ismember(gene_set,column_expr_genes);
G2=G2(:,b(a==1));
column_expr_genes=column_expr_genes(b(a==1));

C_feasible = double(C_adjacency>0);

if pre_synaptic==1
    C_chemical=C_chemical';
    C_feasible=C_feasible';
end

if ~isempty(neuron_2)
    close all;sig_genes=expression_test_neuron_vs_neuron(neuron_1,neuron_2,G1,row_expr_genes,chemical_conn_matrix_row_neurons,double(C_chemical>0),C_chemical,include_self_gene_expr,unique_partners,0,number_of_heterophilic_genes);
else
    close all;sig_genes=expression_test_neuron_partners_vs_adjacent(neuron_1,C_feasible,G1,row_expr_genes,chemical_conn_matrix_row_neurons,double(C_chemical>0),C_chemical,include_self_gene_expr,unique_partners,0,number_of_heterophilic_genes);
end