% 
% Function to parse the position values from the MOPSO into integer values
% represening the actual number of blocks for each peer
%

function Blocks=ParsePositionToBlocks(pop_position_peers, blocks_min, blocks_max)
    
    Blocks_float=blocks_min+(blocks_max-blocks_min).*pop_position_peers;
    
    for pb = 1:numel(pop_position_peers)
        if Blocks_float(pb) > blocks_max 
            Blocks_float(pb) = 200;
        end

        if Blocks_float(pb) < blocks_min 
            Blocks_float(pb) = 1;
        end
    end
    
    Blocks = round(Blocks_float);

end