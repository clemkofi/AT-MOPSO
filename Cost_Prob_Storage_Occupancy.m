%
% Objective Functions for all blocks 
%

function [z, pop_block_peers] = Cost_Prob_Storage_Occupancy(norm_query_probs, pop_position_peers, numPeers, peer_weights, numBlocks)
    
    pop_block_peers = ParsePositionToBlocks(pop_position_peers, 1, 200); 
    
    P_I = ones(numPeers,1);

    % select the query probability sum for all peers
    for peer = 1:numPeers
        num_blocks_peer = pop_block_peers(peer);
        P_I(peer) = norm_query_probs(num_blocks_peer);
    end
        
    % calculate the storage cost
    k = 0.01;
    size_of_block = 1;
    
    sum_of_peer_blocks = 0;
    
    for peer_block = 1:numPeers
        sum_of_peer_blocks = sum_of_peer_blocks + pop_block_peers(peer_block);
    end
    
    storage_cost = k * size_of_block * sum_of_peer_blocks;
    
    
    % calculate the local space occupancy
    local_space_occupancy = 0;
    
    for peer_space = 1:numPeers
        peer_space_occpancy = (exp((numBlocks-pop_block_peers(peer_space))/numBlocks)-1)/(exp(1) - 1);
        peer_space_occpancy_weighted = peer_space_occpancy * peer_weights(peer_space);
        
        % add it to the local space occupancy
        local_space_occupancy = local_space_occupancy + peer_space_occpancy_weighted;
        
    end
    
    
    z = [P_I
       storage_cost
       local_space_occupancy];

end