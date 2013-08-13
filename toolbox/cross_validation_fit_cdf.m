function P_coverage_cv=cross_validation_fit_cdf(dJ_body,dJ_line,conf_level)

n_trials=length(dJ_body)
P_coverage_cv=zeros(n_trials,1);
for i=1:n_trials
  % extract the dJ_body values for all the trials except trial i
  dJ_body_tao_flat=zeros(0,1);
  for j=1:n_trials
    if j~=i
      dJ_body_tao_flat=[dJ_body_tao_flat dJ_body{j}];  %#ok
    end
  end
  
  % fit the params of a scaled chi2 distro
  x_bar=mean(dJ_body_tao_flat);
  s2=var(dJ_body_tao_flat);
  a=x_bar;
  dof=2*x_bar^2/s2;
  cdf_dJ_emp=chi2scaledcdf(dJ_line,dof,a);
  
  % extract the empirical CDF for that
  %cdf_dJ_emp=empirical_cdf(dJ_body_tao_flat,dJ_line);
  
  % entries in cdf_dJ_emp have to be unique for interpolation
  [cdf_dJ_emp_unique,i_keep]=unique(cdf_dJ_emp);
  dJ_line_unique=dJ_line(i_keep);
  dJ_crit=interp1(cdf_dJ_emp_unique,dJ_line_unique,conf_level);
  
  % calculate the empirical coverage prob for this replicate
  P_coverage_cv(i)=mean(dJ_body{i}<dJ_crit);  
end

end
