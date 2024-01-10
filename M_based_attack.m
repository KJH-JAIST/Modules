function Removelist = M_based_attack(G, Module_number)

tempG = G;
Edgelist = tempG.Edges;
Inter_nodes = [];
Removelist = [];
%%
for i = 1:size(Edgelist,1)
    e = Edgelist{i,['EndNodes']};
    u = str2double(e(1));
    v = str2double(e(2));
    if Module_number(u) ~= Module_number(v)
        Inter_nodes = [Inter_nodes, u, v];
    end
end
Inter_nodes = unique(Inter_nodes);
BC_Inter_G = centrality(tempG, 'betweenness');
BC_Inter = BC_Inter_G(Inter_nodes);

tempsort = [Inter_nodes' BC_Inter randperm(length(Inter_nodes))'];
tempsort = sortrows(tempsort, [2,3], 'descend');
attack_seq = tempsort(:,1);

for i = 1:length(attack_seq)
    remove_candidate = attack_seq(i);
    candidate_ID = findnode(tempG, string(remove_candidate));
    [bin, binsize] = conncomp(tempG);
    idx = binsize(bin) == max(binsize);
    SG = subgraph(tempG, idx);
    if findnode(SG, string(remove_candidate)) == 0
        continue
    end

    NN_candidate = neighbors(tempG, string(remove_candidate));
    isInter = 0;
    Mn_candidate = Module_number(candidate_ID);
    for j = 1:length(NN_candidate)
        if Mn_candidate ~= Module_number(str2double(NN_candidate(j)))
            isInter = 1;
        end
    end
    if isInter == 0
        continue
    else
        tempG = rmnode(tempG, string(attack_seq(i)));
        Removelist = [Removelist, attack_seq(i)];
    end
end

rejected = 1:1:numnodes(G);
rejected = setdiff(rejected, Removelist);
BC_rejected = BC_Inter_G(rejected);

tempsort = [rejected' BC_rejected randperm(length(rejected))'];
tempsort = sortrows(tempsort, [2,3], 'descend');
rejected_seq = tempsort(:,1);

Removelist = [Removelist, rejected_seq'];

end