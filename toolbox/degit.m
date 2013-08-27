function path_no_git=degit(path_raw)
% eliminate .git directories from a path string

path_raw_as_array=split_path(path_raw);
path_no_git_as_array=cell(0,1);
for i=1:length(path_raw_as_array)
  k=strfind(path_raw_as_array{i},'.git');
  if isempty(k)
    path_no_git_as_array{end+1}=path_raw_as_array{i};  %#ok
  end
end
path_no_git=combine_path(path_no_git_as_array);

end

