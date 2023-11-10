function config_name = set_test_config(data_path, ckpt, data_size)
%SET_TES_CONFIG 此处显示有关此函数的摘要
%   此处显示详细说明
config_name = fullfile(data_path,'test_config.yaml');

fid=fopen(config_name,'wt');

fprintf(fid,'model_path: %s\n', ckpt);

fprintf(fid,'model:\n');
fprintf(fid,'  name: UNet3D\n');
fprintf(fid,'  in_channels: 3\n');
fprintf(fid,'  out_channels: 4\n');
fprintf(fid,'  layer_order: gcr\n');
fprintf(fid,'  f_maps: [32, 64, 128, 256]\n');
fprintf(fid,'  num_groups: 8\n');
fprintf(fid,'  final_sigmoid: true\n');
fprintf(fid,'  is_segmentation: true\n');

fprintf(fid,'predictor:\n');
fprintf(fid,'  name: ''StandardPredictor''\n');

fprintf(fid,'loaders:\n');
fprintf(fid,'  batch_size: 1\n');
fprintf(fid,'  mirror_padding: [0, 0, 0]\n');
fprintf(fid,'  raw_internal_path: raw\n');
fprintf(fid,'  num_workers: 0\n');
fprintf(fid,'  test:\n');
fprintf(fid,'    file_paths:\n');
fprintf(fid,'      - %s\n', data_path);
fprintf(fid,'    slice_builder:\n');
fprintf(fid,'      name: SliceBuilder\n');
fprintf(fid,'      patch_shape: [%d, %d, %d]\n',data_size(3),data_size(2),data_size(1));
fprintf(fid,'      stride_shape: [100, 160, 160]\n');
fprintf(fid,'    transformer:\n');
fprintf(fid,'        raw:\n');
fprintf(fid,'          - name: Standardize\n');
fprintf(fid,'          - name: ToTensor\n');
fprintf(fid,'            expand_dims: true\n');

fclose(fid);
end

