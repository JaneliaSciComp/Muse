function [syl_name,i_start,i_end,f_lo,f_hi,r_head_pels,r_tail_pels]= ...
  read_voc_index(voc_index_file_name)

% load the index
s=load(voc_index_file_name);
voc_index=s.mouse;
clear s;

% first things first
n_voc=length(voc_index);

% figure out what form the position data is in
spanish=false;
chinese=false;
if isfield(voc_index,'pos_data')
  % this is how the older files store the position data (pre April 2013)
  if n_voc>=1
    pos=voc_index(1).pos_data;
    n_mice=length(pos);
    if n_mice>0
      spanish=isfield(pos(1),'x_head');
    end
  end
else
  % newer files (post April 2013) store the position data in a separate file
  [voc_index_dir_name,voc_index_base_name,voc_index_extension]=fileparts(voc_index_file_name);
  pos_index_file_name=fullfile(voc_index_dir_name,[voc_index_base_name '_Position' voc_index_extension]);
  s=load(pos_index_file_name);
  pos_index=s.mouse_position;
  chinese=true;
  if n_voc>=1
    pos=pos_index(1);
    n_mice=length(pos.pos_data_nose_x);
  else
    n_mice=0;
  end
end

% prep the arrays
syl_name=cell(n_voc,1);
i_start=zeros(n_voc,1);
i_end=zeros(n_voc,1);
f_lo=zeros(n_voc,1);
f_hi=zeros(n_voc,1);
r_head_pels=zeros(2,n_voc,n_mice);
r_tail_pels=zeros(2,n_voc,n_mice);

% iterate over the vocalizations
for i=1:n_voc
  % read in the syllable name
  syl_name{i}=voc_index(i).syl_name;

  % read the start, end index
  i_start(i)=floor(voc_index(i).start_sample_fine);
  i_end(i)=ceil(voc_index(i).stop_sample_fine);

  % read the band where the vocalization is
  f_lo(i) = voc_index(i).lf_fine;  % Hz
  f_hi(i) = voc_index(i).hf_fine;

  % get head & tail positions, in pels
  if spanish
    pos=voc_index(i).pos_data;
    for j=1:n_mice
      r_head_pels(:,i,j)=[pos(j).x_head ...
                          pos(j).y_head]';  % pels
      r_tail_pels(:,i,j)=[pos(j).x_tail ...
                          pos(j).y_tail]';  % pels
    end
  elseif chinese
    pos=pos_index(i);
    for j=1:n_mice
      r_head_pels(:,i,j)=[pos.pos_data_nose_x(j) ...
                          pos.pos_data_nose_y(j)]';  % pels
      r_tail_pels(:,i,j)=[pos.pos_data_tail_x(j) ...
                          pos.pos_data_tail_y(j)]';  % pels
    end    
  else
    % "english"
    pos=voc_index(i).pos_data;
    for j=1:n_mice
      r_head_pels(:,i,j)=[pos(j).nose_x ...
                          pos(j).nose_y]';  % pels
      r_tail_pels(:,i,j)=[pos(j).tail_x ...
                          pos(j).tail_y]';  % pels
    end
  end
end

% % filter out the "vocalizations" with a low-freq cutoff less than 30 kHz
% % also, get rid of those with f_lo>f_hi, they look like just noise
% keep=(f_lo>30e3)&(f_lo<f_hi);
% i_syl=i_syl(keep);
% i_start=i_start(keep);
% i_end=i_end(keep);
% f_lo=f_lo(keep);
% f_hi=f_hi(keep);
% r_head_pels=r_head_pels(:,keep,:);
% r_tail_pels=r_tail_pels(:,keep,:);

end
