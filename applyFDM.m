% Function to apply EWT
function rhythm_data = applyFDM(data, num_cells, num_channels, Fs, freq_bands, trans_type)
    rhythm_data = cell(num_cells, 1);
    for i = 1:num_cells
        channel_results = cell(num_channels, 1);
        for j = 1:num_channels
            channel_data = data{i, 1}(j, :);
            t1 = (0:length(channel_data) - 1) / Fs;
            reconstructed_IMFs = sp_DFTOrthogonalOrFIR_IIR_LINOEP(channel_data, t1, Fs, freq_bands, trans_type);
            channel_results{j} = reconstructed_IMFs;
        end
        rhythm_data{i} = channel_results;
    end
end