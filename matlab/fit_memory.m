function k = fit_memory(m_short, m_long, t_short, t_long)
% Least squares fit of m(t) = exp(-k t) through m(0) = 1 and two washout points.
  ys = log(max(m_short, 1e-6));
  yl = log(max(m_long, 1e-6));
  k = -(t_short .* ys + t_long .* yl) ./ (t_short.^2 + t_long.^2);
  bad = (m_short <= 0.02) | (m_long <= 0.02) | (m_short > 1.5) | (m_long > 1.5) | (k <= 0);
  k(bad) = NaN;
end
