function [G_attacked,removelist] = HBCA_Non_renewal(G,c)
n = numnodes(G); % Size of graph
d = round(n*c); % the number of removed nodes
L = zeros(2,n);
Names = string(G.Nodes.Name);
L(1,:) = Names';
L(2,:) = centrality(G,'betweenness');
removelist = zeros(1,d);

for i = 1:d
    d_max = max(L(2,:));
    idx = find(L(2,:)==d_max);
    if length(idx) > 1
        idx = datasample(idx,1);
    end
    removelist(i) = string(L(1,idx));
    L(:,idx) = [];
end
G_attacked = rmnode(G,removelist);

end

