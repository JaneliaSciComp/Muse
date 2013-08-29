function [r_head_fake,r_tail_fake]=random_mouse_locations(R,r_head,r_tail,n_fake_mice)

R_min=min(R(1:2,:),[],2);  % drop z coord
R_max=max(R(1:2,:),[],2);
are_all_heads_and_tails_inside_bounding_box=false;
while ~are_all_heads_and_tails_inside_bounding_box ,
  f_rand=rand(2,n_fake_mice);
  r_head_fake=bsxfun(@times,1-f_rand,R_min)+bsxfun(@times,f_rand,R_max);
  length=norm(r_head-r_tail);
  theta_rand=2*pi*rand(1,n_fake_mice);
  theta_rand_hat=[cos(theta_rand);sin(theta_rand)];
  r_tail_fake=bsxfun(@plus,r_head_fake,length*theta_rand_hat);
  is_tail_inside_bounding_box = ...
    all(bsxfun(@lt,R_min,r_tail_fake),1) & ...
    all(bsxfun(@lt,r_tail_fake,R_max),1)   % 1 x n_fake_mice
  are_all_heads_and_tails_inside_bounding_box = ...
    all(is_tail_inside_bounding_box);  % heads are automatically inside
end

end
