% 
% Function to calculate the query propability of all the blocks for a peer 
% === For the fixed query frequency case
%

function [normalized_sum_query_prob, normalized_query_prob] = CalculateQueryProbability_FFixed(F0, numBlocks)

    %% query probability of the blocks
    query_prob = F0.* ones(numBlocks,1);
    normalized_query_prob = zeros(numBlocks,1);
    normalized_sum_query_prob = zeros(numBlocks,1);
   
    %% function for the query frequency used in the integral
    prob_freq_fun = @(t) (F0.*t.^0);
    
    sum_query_prob = 0;
    
    %% computing the query probabilities for each block
    for j = 1:numBlocks
        if j == numBlocks
            query_prob(j) = F0;
        else
            
            % find the integral of the query frequency for a particular block number j 
            prob_freq = integral(prob_freq_fun,0,numBlocks-j);
            
            % calculate the query probability a particular block number j
            query_prob(j) = 1/prob_freq; 
        
        end
        
        % add the computed query probability to the sum of query
        % probabilities
        sum_query_prob = sum_query_prob + query_prob(j);
    end
    
    %% computing the normalized query probabilities for each block
    for k = 1:numBlocks
        if k == numBlocks
            normalized_query_prob(k) = F0/sum_query_prob;
        else
            
            % find the integral of the query frequency for a particular block number j 
            norm_prob_freq = integral(prob_freq_fun,0,numBlocks-k);
            
            % calculate the query probability a particular block number j
            normalized_query_prob(k) = 1/(sum_query_prob*norm_prob_freq); 
        
        end
        
        % compute the sum at each block number and store it
        if k == 1
            normalized_sum_query_prob(k) = normalized_query_prob(k);
        else
            normalized_sum_query_prob(k) = normalized_query_prob(k) + normalized_sum_query_prob(k-1);
        end
       
    end
    
end