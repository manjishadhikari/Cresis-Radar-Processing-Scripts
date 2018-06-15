%Run save data features
if ~exist('gRadar','var')
  startup;
end

param_override = [];
param_override.sched.type='custom_torque';
param_override.sched.type = 'no scheduler';
%param_override.sched.cluster_size = Inf;
param_override.sched.rerun_only = false;
param_override.sched.submit_arguments    = '-l nodes=1:ppn=1,walltime=00:60:00';
%param_override.sched.stop_on_fail = true;
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%=======ENTER LINE NUMBER TO PROCESS=============
param_override.location={'Peterman'};    %Peterman or Jacobshavn
      %Set to 1 to save figure only; 0 to save data features
param_override.lines=[1:8];    %eg {1:36 jacob vl   1:25 cl) { Peterman==> {1:29} crossline and {1:20}vertical}

%==================================================

calculate_reflectivity_tracker(param_override)
return