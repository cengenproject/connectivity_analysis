function [sig_genes,beta,pval,gene_cross]=expression_test_neuron_vs_neuron(neuron_1,neuron_2,G,gene_set,neuron_set,connectome,weighted_connectome,withself,unique,normalization,numhetero)
if nargin<11
    withself=1;
    unique=0;
    normalization=1;
    numhetero=20;
end

numhetero=min(numhetero,length(gene_set));
set(groot, 'DefaultTextInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','none');
i=find(strcmpi(neuron_set,neuron_1));
j=find(strcmpi(neuron_set,neuron_2));

partners_1 = find(connectome(i,:));
partners_2 = find(connectome(j,:));
partners_common=intersect(partners_1,partners_2);

partners_1=[setdiff(partners_1,partners_common) partners_common];
partners_2=[fliplr(partners_common) setdiff(partners_2,partners_common)];

if unique==1
partners_1_tmp = setdiff(partners_1,partners_2);
partners_2_tmp = setdiff(partners_2,partners_1);

partners_1=partners_1_tmp;
partners_2=partners_2_tmp;
end
for u=1:length(gene_set)
    for v=1:length(gene_set)
        if withself==1
            gene_interaction_1=log1p(G(partners_1,u)*G(i,v));
            gene_interaction_2=log1p(G(partners_2,u)*G(j,v));
        else
            gene_interaction_1=log1p(G(partners_1,u));
            gene_interaction_2=log1p(G(partners_2,u));
        end
        if normalization==1
            gene_interaction_1 = -log1p(weighted_connectome(partners_1,i)) + gene_interaction_1 ;
            gene_interaction_2 = -log1p(weighted_connectome(partners_2,j)) + gene_interaction_2 ;
        end
        [h,p,ci,stats]=ttest2(gene_interaction_1,gene_interaction_2);
        logfold=mean(gene_interaction_1,1)-mean(gene_interaction_2,1);
        beta(u,v)=logfold;
        tstat(u,v)=stats.tstat;
        pval(u,v)=p;
        [u v]
    end
end

% [~,idx]=sort(diag((beta)),'ascend','MissingPlacement','last');

[~,idx_low]=sort(diag((beta)),'ascend','MissingPlacement','last');
[~,idx_high]=sort(diag((beta)),'descend','MissingPlacement','last');

idx=[idx_low(1:numhetero)' flipud(idx_high(1:numhetero))'];

figure('units','normalized','outerposition',[0 0 0.5 1])


if withself==1
    if normalization==1
        imagesc([log1p(G(partners_1,idx).*G(i,idx))-log1p(weighted_connectome(partners_1,i));log1p(G(partners_2,idx).*G(j,idx))-log1p(weighted_connectome(partners_2,j))]);colormap(othercolor('Blues9'));
    else
    imagesc([log1p(G(partners_1,idx).*G(i,idx));log1p(G(partners_2,idx).*G(j,idx))]);colormap(othercolor('Blues9'));
    end
else
    if normalization==1
        imagesc([log1p(G(partners_1,idx))-log1p(weighted_connectome(partners_1,i));log1p(G(partners_2,idx))-log1p(weighted_connectome(partners_2,j))]);colormap(othercolor('Blues9'));
    else
    imagesc([log1p(G(partners_1,idx));log1p(G(partners_2,idx))]);colormap(othercolor('Blues9'));
    end
end
colorbar
a=neuron_set(partners_1)';
for t=1:length(a)
    a{t}=[neuron_1 '+' a{t}];
end
hold on
if and(~isempty(partners_common),unique==0)
plot([1 2*numhetero],[length(partners_1)-length(partners_common)+0.5 length(partners_1)-length(partners_common)+0.5],'m--','LineWidth',2);
plot([1 2*numhetero],[length(partners_1)+length(partners_common)+0.5 length(partners_1)+length(partners_common)+0.5],'m--','LineWidth',2);
end
plot([1 2*numhetero],[length(a)+0.5 length(a)+0.5],'r-','LineWidth',2);
b=neuron_set(partners_2)';
for t=1:length(b)
    b{t}=[neuron_2 '+' b{t}];
end
a=a';
b=b';
plot([numhetero+0.5 numhetero+0.5],[1 length(partners_1)+length(partners_2)],'r-','LineWidth',2);
set(gca,'Ytick',(1:size([gene_interaction_1;gene_interaction_2],1)),'YTickLabel',[a;b],'Xtick',(1:length(gene_set(idx))),'XTickLabel',gene_set(idx),'XTickLabelRotation',90,'FontWeight','bold');
title('Homophilic CAMs -- above red line: neuron-1+synaptic partners, below red line: neuron-2 + synaptic partners');


% figure('units','normalized','outerposition',[0.5 0 0.5 1])
% 
% dB=diag(beta);
% dP=diag(pval);
% plot(dB,-log10(dP),'k.','MarkerSize',10);
% hold on
% idx=find(and(dB<0,dP<=0.05/length(gene_set)));
% plot(dB(idx),-log10(dP(idx)),'r.','MarkerSize',14);
% text(dB(idx),-log10(dP(idx)),gene_set(idx),'Color','k','FontSize',14,'FontWeight','bold');
% idx=find(and(dB>0,dP<=0.05/length(gene_set)));
% plot(dB(idx),-log10(dP(idx)),'b.','MarkerSize',14);
% text(dB(idx),-log10(dP(idx)),gene_set(idx),'Color','k','FontSize',14,'FontWeight','bold');
% xlabel('Log-fold change')
% ylabel('-log10 p-val');
% axis square
% title('Homophilic CAMs volcano plot')
% legend({'',['Significantly more in ' neuron_set{j} ' and its neighbors'],['Significantly more in ' neuron_set{i} ' and its neighbors']});
% 
% sig_genes=gene_set(dP<=0.05/length(gene_set));

figure('units','normalized','outerposition',[0.5 0 0.5 1])

dB=beta(:);
dP=pval(:);
[A,B]=meshgrid(gene_set,gene_set);
gene_set_cross=cellstr([char(A(:)) repmat('x ',size(A(:))) char(B(:))]);
plot(dB,-log10(dP),'k.','MarkerSize',20);
hold on
idx=find(and(dB<-0.2,dP<=0.05/(length(gene_set_cross))));
plot(dB(idx),-log10(dP(idx)),'r.','MarkerSize',20);
text(dB(idx),-log10(dP(idx)),gene_set_cross(idx),'Color','k','FontSize',20,'FontWeight','bold');
idx=find(and(dB>0.2,dP<=0.05/(length(gene_set_cross))));
plot(dB(idx),-log10(dP(idx)),'b.','MarkerSize',20);
text(dB(idx),-log10(dP(idx)),gene_set_cross(idx),'Color','k','FontSize',20,'FontWeight','bold');
xlabel('Log-fold change')
ylabel('-log10 p-val');
axis square
title('Homophilic + heterophilic CAMs volcano plot')
legend({'',['Significantly more in ' neuron_set{j} ' and its neighbors'],['Significantly more in ' neuron_set{i} ' and its neighbors']});
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
sig_genes=gene_set_cross(dP<=0.05/length(gene_set_cross));

% [~,idx_low]=sort(-log10(dP).*(sign(-dB)),'descend','MissingPlacement','last');
[~,idx_low]=sort(dB,'ascend','MissingPlacement','last');
[subi_low,subj_low]=ind2sub(size(pval),idx_low);

% [~,idx_high]=sort(-log10(dP).*(sign(dB)),'descend','MissingPlacement','last');
[~,idx_high]=sort(dB,'descend','MissingPlacement','last');
[subi_high,subj_high]=ind2sub(size(pval),idx_high);


idx_a=[subi_low(1:numhetero)' flipud(subi_high(1:numhetero))'];
idx_b=[subj_low(1:numhetero)' flipud(subj_high(1:numhetero))'];

gene_set_hetero=gene_set_cross([idx_low(1:numhetero)' flipud(idx_high(1:numhetero))']);

figure('units','normalized','outerposition',[0.5 0 0.5 1])
if withself==1
    if normalization==1
        imagesc([log1p(G(partners_1,idx_a).*G(i,idx_b))-log1p(weighted_connectome(partners_1,i));log1p(G(partners_2,idx_a).*G(j,idx_b))-log1p(weighted_connectome(partners_2,j))]);colormap(othercolor('Blues9'));
    else
    imagesc([log1p(G(partners_1,idx_a).*G(i,idx_b));log1p(G(partners_2,idx_a).*G(j,idx_b))]);colormap(othercolor('Blues9'));
    end
else
    if normalization==1
        imagesc([log1p(G(partners_1,idx_a))-log1p(weighted_connectome(partners_1,i));log1p(G(partners_2,idx_a))-log1p(weighted_connectome(partners_2,j))]);colormap(othercolor('Blues9'));
    else
    imagesc([log1p(G(partners_1,idx_a));log1p(G(partners_2,idx_a))]);colormap(othercolor('Blues9'));
    end
end
colorbar
a=neuron_set(partners_1)';
for t=1:length(a)
    a{t}=[neuron_1 '+' a{t}];
end
hold on
if and(~isempty(partners_common),unique==0)
plot([1 2*numhetero],[length(partners_1)-length(partners_common)+0.5 length(partners_1)-length(partners_common)+0.5],'m--','LineWidth',2);
plot([1 2*numhetero],[length(partners_1)+length(partners_common)+0.5 length(partners_1)+length(partners_common)+0.5],'m--','LineWidth',2);
end
plot([1 length(gene_set_hetero)],[length(a)+0.5 length(a)+0.5],'r-','LineWidth',2);
plot([numhetero+0.5 numhetero+0.5],[1 length(partners_1)+length(partners_2)],'r-','LineWidth',2);
b=neuron_set(partners_2)';
for t=1:length(b)
    b{t}=[neuron_2 '+' b{t}];
end
a=a';
b=b';
set(gca,'Ytick',(1:size([gene_interaction_1;gene_interaction_2],1)),'YTickLabel',[a;b],'Xtick',(1:length(gene_set_hetero)),'XTickLabel',gene_set_hetero,'XTickLabelRotation',90,'FontWeight','bold');
title('Top 20 Heterophilic CAMs -- above red line: neuron-1+synaptic partners, below red line: neuron-2 + synaptic partners');
