

% Author: Hanieh Soleimani


% function my_pca: computes the PCA of two given dataset:
% Input: dataset, the original data set.
%        v, the number of eigenvectors.
% Output: reduction, reduced data from the original dataset 
%         err, is the error of the PCA process
%         reconstructed_data 

function [reduction, err, reconstructed_data] = my_pca(file, v)

% load the data set given in the assignmnet
data_set = importdata(file);


% - subtracting mean value from the data set
% - compute the covariance matrix
% - finding the eigenvalues and corresponding eigenvectors.
mean_value = sum(data_set)/ size(data_set,1);
dataSubMean = bsxfun(@minus, data_set , mean_value);
cov_matrix = (dataSubMean' * dataSubMean)/ (size(data_set,1)- 1);
[ V, D] = eig(cov_matrix);

%%%%%%%%%%%%%%%%%%%%%%% PLOT THE EIGENVALUES %%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
semilogy(size(data_set,2):-1:1,diag(D))
xlim([0,10]);
title('Scree Plot');
xlabel('Component Number');
ylabel('Eigenvalue');



%%%%%%%%%%%%%%%%%%%%%%%   PC COMPONENT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  first step:  sorting the eigenvalues in descending order,
%  second step: extract the eigenvectors corresponded to the sorted eigenvalues
[~, sorted_D] = sort(diag(D), 'descend'); 
sorted_V = V(:,sorted_D);


% finding the projection matrix
proj_matrix = sorted_V(:,1:v);
reduction = data_set * proj_matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%   PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
histogram(data_set);
hold on;

% reconstructing the original data set and plot the comparison plot
reconstructed_data = (reduction * proj_matrix') + bsxfun(@plus, data_set , mean_value);
histogram(reconstructed_data);


title('Comparing Original with Reconstructed data');
xlabel('values');
ylabel('Frequency');


figure(3)
if(v==3)
    scatter3(reduction(:,1),reduction(:,2), reduction(:,3));
else 
     if (v==2)
        scatter(reduction(:,1),reduction(:,2));
    end
end

% plot of reduced data into 2 dimensions
title('reduced\_data');
xlabel(' First\_Component');
ylabel('second\_Component');
% plot for the 3rd component
if(v==3)
    zlabel('Third\_Component');
end



%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% using the formula on the note on the icorsi %%%%%%%%%%%%%%%%

d_v = diag(D);
non_used_eigenval = d_v(v+1:size(data_set,1));
numerator = sum(non_used_eigenval,2);

denominator = sum(d_v(1:size(data_set,1)));

err = (1 - (numerator/denominator))*100;




end 






