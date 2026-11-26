% Independent Octave/MATLAB implementation used as a cross-language check in CI.
% Usage: octave --no-gui --eval "run_all('build/octave')"
function run_all(folder)
  addpath(fileparts(mfilename('fullpath')));
  p_on = struct('r_S', -0.8, 'r_P', 0.0, 'a', 0.05, 'b', 0.0);
  p_off = struct('r_S', 0.6, 'r_P', 0.1, 'a', 0.001, 'b', 0.3);
  traj = switching_model(p_on, p_off, 21, 21, 0.25);
  M = [p_on.r_S - p_on.a, p_on.b; p_on.a, p_on.r_P - p_on.b];
  lam_on = max(real(eig(M)));
  M = [p_off.r_S - p_off.a, p_off.b; p_off.a, p_off.r_P - p_off.b];
  lam_off = max(real(eig(M)));
  fid = fopen(fullfile(folder, 'model_octave.csv'), 'w');
  fprintf(fid, 'quantity,value\n');
  fprintf(fid, 'lambda_on,%.10g\nlambda_off,%.10g\n', lam_on, lam_off);
  fprintf(fid, 'total_day21,%.10g\ntotal_day42,%.10g\n', sum(traj(85, 2:3)), sum(traj(end, 2:3)));
  fprintf(fid, 'persister_fraction_day21,%.10g\n', traj(85, 3) / sum(traj(85, 2:3)));
  fclose(fid);
  f = fullfile(folder, 'memory_input.csv');
  if exist(f, 'file')
    fid = fopen(f, 'r');
    header = fgetl(fid);
    C = textscan(fid, '%s %s %f %f %f', 'Delimiter', ',');
    fclose(fid);
    k = fit_memory(C{3}, C{4}, 1.0, C{5});
    fid = fopen(fullfile(folder, 'memory_octave.csv'), 'w');
    fprintf(fid, 'line,gene_id,k\n');
    for i = 1:numel(k)
      fprintf(fid, '%s,%s,%.10g\n', C{1}{i}, C{2}{i}, k(i));
    end
    fclose(fid);
  end
  disp('octave cross-check written');
end
