function [pid, ppid] = print_pointwise_pid_iccs(Pxxy,normalise)

if nargin<2
    normalise = false;
end

[Na Nb Ny] = size(Pxxy);
P2 = marg_maxent2(Pxxy);

Pxx = sum(Pxxy,3);
Pxx2 = sum(P2,3);
Py = squeeze(sum(sum(Pxxy,1),2));
Pay = squeeze(sum(Pxxy,2));
Pby = squeeze(sum(Pxxy,1));
Pa = squeeze(sum(sum(Pxxy,3),2));
Pb = squeeze(sum(sum(Pxxy,3),1));

fprintf(1,'\n')
% fprintf(1,'     Pxxy      : \n\n')
fprintf(1,' Paby   A B Y  : i(a;y) i(b;y) i(ab;y) :   red    unqA   unqB   syn \n')

pid = zeros(1,4);
ppid = zeros([size(Pxxy) 4]);
for ai=1:Na
    for bi=1:Nb
        for yi=1:Ny
            if Pxxy(ai,bi,yi)==0
                continue
            end
             
            dsj = log2( Pxxy(ai,bi,yi) / (Pxx(ai,bi)*Py(yi)) );
            dsj2 = log2( P2(ai,bi,yi) / (Pxx2(ai,bi)*Py(yi)) );
            dsa = log2( Pay(ai,yi) ./ (Pa(ai)*Py(yi)) );
            dsb = log2( Pby(bi,yi) ./ (Pb(bi)*Py(yi)) );
            
            % NB this is new not in Iccs
            % because in some terms the P2 maxent removes what looks like
            % it should be redundnat misinformation idea was to tkae the
            % max (unsigned) overlap of both P and P2
            % this lets P2 induce sharing when there is
            overlapP2 = dsa + dsb - dsj2;       
            if sign(dsa)~=sign(overlapP2) || sign(dsa)~=sign(dsb) || sign(dsa)~=sign(dsj2)
                overlapP2 = 0;
            end
            overlapPorig = dsa + dsb - dsj;
            if sign(dsa)~=sign(overlapPorig) || sign(dsa)~=sign(dsb) || sign(dsa)~=sign(dsj)
                overlapPorig = 0;
            end
            
            overlaps = [overlapP2 overlapPorig];
            [~,idx] = max(abs(overlaps));
            overlap = overlaps(idx);
            unqA = dsa - overlap;
            unqB = dsb - overlap;
            syn = dsj - unqA - unqB - overlap;
            
            pid(1) = pid(1) + Pxxy(ai,bi,yi) * overlap;
            pid(2) = pid(2) + Pxxy(ai,bi,yi) * unqA;
            pid(3) = pid(3) + Pxxy(ai,bi,yi) * unqB;
            pid(4) = pid(4) + Pxxy(ai,bi,yi) * syn;
            % could use P2 for overlap, but then sum up accroding to
            % original dsitribution
            if normalise
                overlap = Pxxy(ai,bi,yi) * overlap;
                unqA = Pxxy(ai,bi,yi) * unqA;
                unqB = Pxxy(ai,bi,yi) * unqB;
                syn = Pxxy(ai,bi,yi) * syn;
                dsa = Pxxy(ai,bi,yi) * dsa;
                dsb = Pxxy(ai,bi,yi) * dsb;
                dsj = Pxxy(ai,bi,yi) * dsj;
            end

            ppid(ai,bi,yi,1) = overlap;
            ppid(ai,bi,yi,2) = unqA;
            ppid(ai,bi,yi,3) = unqB;
            ppid(ai,bi,yi,4) = syn;
            fprintf(1,'%6.3f [%d %d %d] : %6.3f %6.3f %6.3f  : %6.3f %6.3f %6.3f %6.3f \n',Pxxy(ai,bi,yi),ai-1,bi-1,yi-1,dsa,dsb,dsj,overlap,unqA,unqB,syn) 
            
        end
    end
end
fprintf('\n')