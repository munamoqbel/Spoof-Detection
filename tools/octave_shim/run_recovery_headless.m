pkg load statistics
cd(getenv('SHIM')); addpath('/home/user/Spoof-Detection'); warning('off', 'all');
run_recovery;
viol = find(err_rec > max(3*sig_rss, 1e-9));
fprintf('violating epochs: %s | state there: %s | err=%.4f 3sig=%.4f\n', mat2str(viol), mat2str(out.state(viol)), err_rec(viol), 3*sig_rss(viol));
fprintf('t_detect=%s t_anchor=%s t_prob_start=%s t_prob_fail=%s t_handback=%s\n', mat2str(out.ev.t_detect), mat2str(out.ev.t_anchor), mat2str(out.ev.t_prob_start), mat2str(out.ev.t_prob_fail), mat2str(out.ev.t_handback));
