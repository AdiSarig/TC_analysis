function [faces_vis_stat, houses_vis_stat] = run_decoding_one_paradigm(paradigm, data,sub_inds)

for ind = 1:length(sub_inds)
    fprintf('%s: Preparing sub %d out of %d...\n', paradigm, ind, length(sub_inds));
    switch paradigm
        case 'masking'
            % prepare data
            faces_vis{1, ind} = data(sub_inds(ind)).V_f;
            faces_vis{2, ind} = data(sub_inds(ind)).IV_f;
            
            houses_vis{1, ind} = data(sub_inds(ind)).V_h;
            houses_vis{2, ind} = data(sub_inds(ind)).IV_h;

%             blank_vis{1, ind} = data.V_b(sub_inds(ind)).EEG;
%             blank_vis{2, ind} = data.IV_b(sub_inds(ind)).EEG;
        case 'ib'
            faces_vis{1, ind} = data.P2_f(sub_inds(ind)).EEG;
            faces_vis{2, ind} = data.P1_f(sub_inds(ind)).EEG;
            
            houses_vis{1, ind} = data.P2_h(sub_inds(ind)).EEG;
            houses_vis{2, ind} = data.P1_h(sub_inds(ind)).EEG;
        otherwise
            % prepare data
            faces_vis{1, ind} = data.V_f(sub_inds(ind)).EEG;
            faces_vis{2, ind} = data.IV_f(sub_inds(ind)).EEG;
            
            houses_vis{1, ind} = data.V_h(sub_inds(ind)).EEG;
            houses_vis{2, ind} = data.IV_h(sub_inds(ind)).EEG;
    end
end

% blank_vis_stat = perm_decoding(blank_vis);
faces_vis_stat = perm_decoding(faces_vis);
houses_vis_stat = perm_decoding(houses_vis);

end

