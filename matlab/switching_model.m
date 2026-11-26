function out = switching_model(p_on, p_off, on_days, off_days, dt)
% Two-state sensitive/persister switching model, matrix exponential stepping.
% p_on, p_off: structs with fields r_S, r_P, a, b (per day).
  if nargin < 5, dt = 0.25; end
  M_on = model_matrix(p_on);
  M_off = model_matrix(p_off);
  n_on = round(on_days / dt);
  n_off = round(off_days / dt);
  x = [0.999; 0.001];
  out = zeros(n_on + n_off + 1, 3);
  out(1, :) = [0, x'];
  E_on = expm(M_on * dt);
  E_off = expm(M_off * dt);
  for i = 1:(n_on + n_off)
    if i <= n_on
      x = E_on * x;
    else
      x = E_off * x;
    end
    out(i + 1, :) = [i * dt, x'];
  end
end

function M = model_matrix(p)
  M = [p.r_S - p.a, p.b; p.a, p.r_P - p.b];
end
