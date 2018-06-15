

function calculate_reflectivity_tracker(param)

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: (%s)\n', dbstack_info(1).name, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');


task_param=param;

% =====================================================================
% Setup the scheduler
% =====================================================================

if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('save_data_features_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
  
elseif ~strcmpi(param.sched.type,'no scheduler')
  fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
    get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
  fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
    get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];
  
  fd = merge_filelists(fd, fd_override);
  
  % Remove SVN files from path
  non_svn_path = zeros(size(fd));
  for fn_idx = 1:length(fd)
    if isempty(strfind(fd{fn_idx},'.svn'))
      non_svn_path(fn_idx) = 1;
    end
  end
  fd = fd(non_svn_path);
  
  % Initialize submission ctrl structure
  global ctrl; % Make this global for convenience in debugging
  ctrl = [];
  ctrl.cmd = 'init';
  ctrl.sched = param.sched;
  ctrl.fd = fd;
  ctrl = create_task(ctrl);
  
  % Prepare submission ctrl structure for queing jobs
  ctrl.cmd = 'task';
end


%%
if 1
for line_idx = 1:length(param.lines)
  task_param.proc_line= param.lines(line_idx);
  
    % Execute tasks/jobs
    fh = @calculate_reflectivity_task;
   
    arg{1} = task_param;
   
    
    if strcmp(param.sched.type,'custom_torque')
     
      create_task_param.conforming = true;
      ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
       fprintf(' sub %d/%d:  Matlab job id: %d. Torque Job id: %d. Processing Line number %d \n',ctrl.batch_id,ctrl.job_id_list(line_idx),line_idx,ctrl.job_id_list(line_idx), task_param.proc_line);
      
    elseif ~strcmp(param.sched.type,'no scheduler')
      [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
     fprintf('Processing line no %d\n',task_param.proc_line);
      retry_fields{job_id,task_id}.frm =  task_param.proc_line;
 
      retry_fields{job_id,task_id}.arg = arg;
   
    else
    fprintf('Processing Line number %d \n', task_param.proc_line);
    [success] = fh(arg{1});
    end
    
  
end

% =======================================================================
% Wait for jobs to complete if a scheduler was used
% =======================================================================
if strcmpi(param.sched.type,'custom_torque')
  % Wait until all submitted jobs to complete
  ctrl = torque_rerun(ctrl);
  if ~all(ctrl.error_mask == 0)
    if ctrl.sched.stop_on_fail
      torque_cleanup(ctrl);
      error('Not all jobs completed, but out of retries (%s)', datestr(now));
    else
      warning('Not all jobs completed, but out of retries (%s)', datestr(now));
      keyboard;
    end
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
   torque_cleanup(ctrl);
  
elseif ~strcmpi(param.sched.type,'no scheduler')
  % ======================================================================
  % Wait for jobs to finish and clean up
  ctrl.cmd = 'done';
  ctrl = create_task(ctrl);
  if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
    % Quit if a bad error occurred
    fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
    if strcmp(ctrl.sched.type,'torque')
      fprintf('Often on the Torque scheduler, these are not bad errors\n');
      fprintf('because of system instabilities (e.g. file IO failure)\n');
      fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
      fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
      keyboard
      if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
        return;
      end
    else
      return
    end
  end
  
  % ======================================================================
  % Retry jobs that failed
  retry = 1;
  while ctrl.error_mask == 2 && retry <= param.sched.max_retries
    fprintf('Tasks failed, retry %d of max %d\n', retry, param.sched.max_retries);
    
    % Bookkeeping (move previous run information to "old_" variables)
    old_ctrl = ctrl;
    old_retry_fields = retry_fields;
    retry_fields = {};
    old_out_recs = out_recs;
    out_recs = {};
    
    % Initialize submission ctrl structure
    ctrl = [];
    ctrl.cmd = 'init';
    ctrl.sched = param.sched;
    ctrl.fd = fd;
    ctrl = create_task(ctrl);
    
    % Prepare submission ctrl structure for queing jobs
    ctrl.cmd = 'task';
    
    % Submit failed tasks, but keep track of these in case they fail again
    for job_idx = 1:length(old_ctrl.jobs)
      for task_idx = old_ctrl.jobs{job_idx}.error_idxs
        [ctrl,job_id,task_id] = create_task(ctrl,fh,2,old_retry_fields{job_idx,task_idx}.arg);
        out_idx = old_retry_fields{job_idx,task_idx}.out_idx;
        fprintf('  %d/%d: Processing records %d to %d in job,task %d,%d (%s)\n', ...
          old_retry_fields{job_idx,task_idx}.frm, old_retry_fields{job_idx,task_idx}.break_idx, ...
          old_out_recs{out_idx}(1), old_out_recs{out_idx}(end), ...
          job_id, task_id, datestr(now));
        retry_fields{job_id,task_id} = old_retry_fields{job_idx,task_idx};
        out_recs{end + 1} = old_out_recs{out_idx};
        retry_fields{job_id,task_id}.out_idx = length(out_recs);
      end
    end
    
    % Wait for tasks to complete and then cleanup
    ctrl.cmd = 'done';
    ctrl = create_task(ctrl);
    retry = retry + 1;
    
    if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
      % Quit if a bad error occurred
      fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
      if strcmp(ctrl.sched.type,'torque')
        fprintf('Often on the Torque scheduler, these are not bad errors\n');
        fprintf('because of system instabilities (e.g. file IO failure)\n');
        fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
        fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
        keyboard
        if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
          return;
        end
      else
        return
      end
    end
  end
  if ctrl.error_mask ~= 0
    fprintf('Not all jobs completed, but out of retries (%s)\n', datestr(now));
    return;
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
end
end
return

