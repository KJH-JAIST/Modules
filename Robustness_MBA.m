function R = Robustness_MBA(G, Removelist)
n = numnodes(G);
Sq = 0;
for i = 1:length(Removelist)
    tempG = G;
    tempG = rmnode(tempG, string(Removelist(1:i)));
    [~, binsize] = conncomp(tempG);
    if isempty(binsize) == 1
    	LCC = 0;
    else
        LCC = max(binsize);
    end
    Sq = Sq + LCC/n;   
end
R = Sq/(n-1);

end