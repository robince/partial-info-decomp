function display_nonzero_terms(lat)

terms = lat.labels(lat.PI>1e-15);
str = sprintf('%s  ',terms{:});
fprintf(1,'NZ PI Terms: %s\n',str);
% fprintf(1,'\n')