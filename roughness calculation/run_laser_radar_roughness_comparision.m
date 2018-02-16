
%Run save data features


param_override = [];
param_override.sched.type='custom_torque';
%param_override.sched.type = 'no scheduler';
%param_override.sched.cluster_size = Inf;
param_override.sched.rerun_only = false;
%param_override.sched.interactive=true;
param_override.sched.submit_arguments    = '-l nodes=1:ppn=1,walltime=00:30:00';
%param_override.sched.stop_on_fail = true;
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%=======ENTER LINE NUMBER TO PROCESS=============
%param_override.location={'Jacobshavn'};    %Peterman or Jacobshavn
param_override.cross_lines_en=1;  %1 for cross line 1 for verticallines 0
param_override.lines=[1];    %eg[1:15]
%param_override.save_en=1;  % 1: Save final results 0:run without saving
%==================================================

laser_radar_roughness_comparision_tracker(param_override)
return