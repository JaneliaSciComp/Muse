function [i_first_tf_rect_in_segment,i_last_tf_rect_in_segment]= ...
  determine_segments(tf_rect_names)

% Find the indices of the first and last t-f rect in each segment.

if isempty(tf_rect_names) ,
  i_first_tf_rect_in_segment=zeros(0,1);
  i_last_tf_rect_in_segment=zeros(0,1);
else
  tf_rect_name_first=tf_rect_names{1};
  if length(tf_rect_name_first)~=17 ,
    error('The "syllable" "names" appear to be in the wrong format.');
  end
  i_segment_for_tf_rect_as_string=cellfun(@(s)s(4:9),tf_rect_names,'UniformOutput',false);
  i_segment_for_tf_rect=str2double(i_segment_for_tf_rect_as_string);
  i_first_tf_rect_in_segment=find(diff([0;i_segment_for_tf_rect])>0);
  i_last_tf_rect_in_segment=[i_first_tf_rect_in_segment(2:end)-1 ; ...
                             length(i_segment_for_tf_rect)];
end

end
