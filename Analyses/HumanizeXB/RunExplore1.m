% RunExplore1.m -- mechanistic sweep: which rate moves which observable (results/explore1.*).
cyc = {'ka','kd','k1','k_1','k2','k_2','kah','kamh'};
srx = {'ksr0','kmsrd','ksrd','ksr2srd','ksrd2sr'};
mk = @(f, v) cell2struct(num2cell(repmat(v, 1, numel(f))), f, 2);
merge = @(a, b) cell2struct([struct2cell(a); struct2cell(b)], [fieldnames(a); fieldnames(b)], 1);
SETS = { 'all_x0.25',        mk(cyc, 0.25)
         'all+srx_x0.25',    merge(mk(cyc, 0.25), mk(srx, 0.25))
         'all+vreg_x0.25',   merge(mk(cyc, 0.25), mk({'v_ref_reg'}, 0.25))
         'vreg_x0.25_only',  mk({'v_ref_reg'}, 0.25)
         'k2_x0.1_rest0.25', merge(mk(setdiff(cyc, {'k2','k_2'}), 0.25), mk({'k2','k_2'}, 0.1))
         'ka_x0.1_rest0.25', merge(mk(setdiff(cyc, {'ka'}), 0.25), mk({'ka'}, 0.1))
         'kah_x0.1_rest0.25', merge(mk(setdiff(cyc, {'kah','kamh'}), 0.25), mk({'kah','kamh'}, 0.1)) };
TAG = 'explore1';
run('C:\home\git\ATP-depletion-and-heart-failure\Analyses\HumanizeXB\RunHumanize.m');
