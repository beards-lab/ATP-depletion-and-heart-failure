fileID = fopen('dymparams.txt','w');
params = getParams([],g);
s = ' =  %3.3f';
d =  datestr(now, 'dd/mm/yy-HH:MM');
fprintf(fileID,['model XBATP_ParamOpt "Optimized parameters, gen at ' d '" \n  extends XBATP(\n']);


fprintf(fileID,['ka ' s ',\n'],params.ka );
fprintf(fileID,['kd ' s ',\n'],params.kd );
fprintf(fileID,['k1 ' s ',\n'],params.k1 );
fprintf(fileID,['km1' s ',\n'],params.k_1);
fprintf(fileID,['k2 ' s ',\n'],params.k2 );
fprintf(fileID,['ksr0  ' s ',\n'],params.ksr0  );
fprintf(fileID,['sigma0' s ',\n'],params.sigma0);
fprintf(fileID,['kmsr  ' s ',\n'],params.kmsr  );
fprintf(fileID,['alpha3' s ',\n'],params.alpha3);
fprintf(fileID,['k3 ' s ',\n'],params.k3 );
fprintf(fileID,['K_T1' s ',\n'],params.K_T1);
fprintf(fileID,['dr' s ',\n'],params.dr);
fprintf(fileID,['kstiff1' s ',\n'],params.kstiff1);
fprintf(fileID,['kstiff2' s ',\n'],params.kstiff2);
fprintf(fileID,['K_T3' s ',\n'],params.K_T3);
fprintf(fileID,['alpha1' s ',\n'],params.alpha1);
fprintf(fileID,['alpha2' s ',\n'],params.alpha2);
fprintf(fileID,['Amax' s ',\n'],params.Amax);
% fprintf(fileID,['mu' s ',\n'],params.mu);
% fprintf(fileID,['kSE' s ',\n'],params.kSE);
fprintf(fileID,['k_passive' s ',\n'],params.k_pas);
% fprintf(fileID,['vmax' s ');\n'],params.vmax);

fprintf(fileID,[');\n'],params.k_pas);
fprintf(fileID,['end XBATP_ParamOpt;\n']);

fclose(fileID);