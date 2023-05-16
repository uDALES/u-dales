%% Setting large scale forcings
% The entries are z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad
% I don't currently use these forcings but if needed add extra lines as 
% appropriate. 


ls = zeros(length(zf), 10);
ls(:,1) = zf;

lscale = fopen([fpath 'lscale.inp.' expnr], 'w');
fprintf(lscale, '%-12s\n', '# SDBL flow');
fprintf(lscale, '%-60s\n', '# z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad');
fprintf(lscale, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-15.9f %-12.6f %-12.6f %-12.6f %-17.12f\n', ls');
fclose(lscale);

