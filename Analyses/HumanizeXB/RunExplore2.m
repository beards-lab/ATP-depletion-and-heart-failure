% RunExplore2.m -- literature-ratio human candidates (results/explore2.*); resumable, re-run to append.
cyc = {'ka','kd','k1','k_1','k2','k_2','kah','kamh'};
mk = @(f, v) cell2struct(num2cell(repmat(v, 1, numel(f))), f, 2);
merge = @(a, b) cell2struct([struct2cell(a); struct2cell(b)], [fieldnames(a); fieldnames(b)], 1);
lit = struct('k2',0.1,'k_2',0.1,'kah',0.1,'kamh',0.1,'k1',0.3,'k_1',0.3,'ka',0.1,'kd',0.25);
litA15 = lit; litA15.ka = 0.15;
v0 = lit; v0.k1 = 0.2; v0.k_1 = 0.2;
v0b = lit; v0b.kd = 0.1;
v0c = lit; v0c.kd = 0.1; v0c.k1 = 0.2; v0c.k_1 = 0.2;
SETS = { 'lit_A',            lit
         'lit_A_ka0.15',     litA15
         'lit_B_vreg',       merge(lit, struct('v_ref_reg',0.1))
         'lit_C_vreg_srx',   merge(merge(lit, struct('v_ref_reg',0.1)), mk({'ksr0','kmsrd','ksrd','ksr2srd','ksrd2sr'}, 0.1))
         'v0_k1x0.2',        v0
         'v0b_kd0.1',        v0b
         'v0c_kd0.1_k1x0.2', v0c
         'uni_0.1',          merge(mk(cyc, 0.1), struct('v_ref_reg',0.1))
         'uni_0.05',         merge(mk(cyc, 0.05), struct('v_ref_reg',0.05)) };
TAG = 'explore2';
run('C:\home\git\ATP-depletion-and-heart-failure\Analyses\HumanizeXB\RunHumanize.m');
