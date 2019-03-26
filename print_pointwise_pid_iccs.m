function print_pointwise_terms(Pxxy)

[Na Nb Ny] = size(Pxxy);
P2 = marg_maxent2(Pxxy);

Pxx = sum(Pxxy,3);
Py = squeeze(sum(sum(Pxxy,1),2));
Pay = squeeze(sum(Pxxy,2));
Pby = squeeze(sum(Pxxy,1));
Pa = squeeze(sum(sum(Pxxy,3),2));
Pb = squeeze(sum(sum(Pxxy,3),1));

fprintf(1,'\n')
% fprintf(1,'     Pxxy      : \n\n')
fprintf(1,' Paby   A B Y  : i(a;y) i(b;y) i(ab;y):   red    unqA   unqB   syn \n')
for ai=1:Na
    for bi=1:Nb
        for yi=1:Ny
            if Pxxy(ai,bi,yi)==0
                continue
            end
            
            dsj = log2( Pxxy(ai,bi,yi) / (Pxx(ai,bi)*Py(yi)) );
            dsj2 = log2( P2(ai,bi,yi) / (Pxx(ai,bi)*Py(yi)) );
            dsa = log2( Pay(ai,yi) ./ (Pa(ai)*Py(yi)) );
            dsb = log2( Pby(bi,yi) ./ (Pb(bi)*Py(yi)) );
            
            overlap = dsa + dsb - dsj2;
            
            if sign(dsa)~=sign(overlap) || sign(dsa)~=sign(dsb) || sign(dsa)~=sign(dsj)
                overlap = 0;
            end
            unqA = dsa - overlap;
            unqB = dsb - overlap;
            syn = dsj - unqA - unqB - overlap;
            fprintf(1,'%6.3f [%d %d %d] : %6.3f %6.3f %6.3f : %6.3f %6.3f %6.3f %6.3f \n',Pxxy(ai,bi,yi),ai-1,bi-1,yi-1,dsa,dsb,dsj,overlap,unqA,unqB,syn) 
            
            
        end
    end
end
fprintf('\n')