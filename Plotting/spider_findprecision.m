function precision  = spider_findprecision(P)

minvals = min(P);
precision = ones(size(minvals));
precision(minvals>10) = 0;
precision(minvals<1) = 2;
precision(minvals<0.01) = 3;
precision(minvals<0.001) = 4;
precision(minvals<0.0001) = 5;

end