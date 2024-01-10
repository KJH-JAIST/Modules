%%
n = 2500;
d = 4;
Num_of_modules = 4;
N = n*Num_of_modules;

W_list = [0.01, 0.05, 0.1, 0.3, 0.5, 0.9];

Q = zeros(1, length(W_list));
R = zeros(1, length(W_list));

ITER = 100;

while 1
    WrongCase = 0;
    Total_A = [];
    Module_number = [];
    
    for i = 1:Num_of_modules
        A = createRandRegGraph(n, d);
        %A = createCompleteGraph(n);
        temp = ones(1,n).*i;
        Module_number = [Module_number, temp];
        Total_A = blkdiag(Total_A,A);
    end

    G = graph(Total_A);
    G.Edges.Weight = [];
    N = length(adjacency(G));
    Numbering = string(1:N); % Attach node labels
    G.Nodes.Name = Numbering';
    
	[~, binsize] = conncomp(G);
    
    if length(binsize) ~= Num_of_modules
        WrongCase = 1;        
    end

    if mean(degree(G)) ~= d
        WrongCase = 1;
    end
    
    if WrongCase == 0
        break;
    end
end

for i = 1:Num_of_modules
    targets = 1:1:n;
    targets = targets + n*(i-1);
    selected_idx1 = targets(randperm(length(targets),1));
    nn_idx1 = neighbors(G, selected_idx1);
    selected_idx2 = nn_idx1(randperm(length(nn_idx1),1));
    G = rmedge(G, selected_idx1, selected_idx2);

    if i == 1
        first_node_idx = selected_idx1;
        prev_back_node = selected_idx2;
    elseif i == Num_of_modules
        G = addedge(G, prev_back_node, selected_idx1);
        G = addedge(G, selected_idx2, first_node_idx);
    else
        G = addedge(G, prev_back_node, selected_idx1);
        prev_back_node = selected_idx2;
    end
end  

Edgelist = G.Edges;
Intra_edges = [];
Inter_edges = [];

for i = 1:size(Edgelist,1)
    e = Edgelist{i,['EndNodes']};
    u = str2double(e(1));
    v = str2double(e(2));
    if Module_number(u) == Module_number(v)
        Intra_edges = [Intra_edges, i];
    else
        Inter_edges = [Inter_edges, i];
    end
end

for k = 1:ITER    
    for j = 1:length(W_list)
        tempG = G;
        temp_intra = Intra_edges;
        W = round(W_list(j)*length(Intra_edges));
        changed_links = 0;
        First_selected_idx = randperm(length(temp_intra),1);
        First_selected_edges = temp_intra(First_selected_idx);
        edge1 = Edgelist{First_selected_edges,['EndNodes']};
        edge1_source = str2double(edge1(1));
        edge1_target = str2double(edge1(2));
        temp_intra(First_selected_idx) = [];
        tempG = rmedge(tempG, string(edge1_source), string(edge1_target));
        First_dest = edge1_source;

        over_ITER = 0;
   
        while changed_links < W
            selected_idx = randperm(length(temp_intra),1);
            selected_edges = temp_intra(selected_idx);
            edge2 = Edgelist{selected_edges,['EndNodes']};
            edge2_source = str2double(edge2(1));
            edge2_target = str2double(edge2(2));

            if (W - changed_links) < 3
                if over_ITER > 500
                    break
                end
                if Module_number(First_dest) == Module_number(edge2_source)
                    over_ITER = over_ITER + 1;
                    continue
                end
                if Module_number(edge2_target) == Module_number(Last_dest)
                    over_ITER = over_ITER + 1;
                    continue
                end
                isedge1 = 0;
                isedge2 = 0;
                isedge3 = 0;
                isedge4 = 0;
                if findedge(tempG, string(First_dest), string(edge2_source)) > 0
                    isedge1 = 1;
                end
                if findedge(tempG, string(First_dest), string(edge2_target)) > 0
                    isedge2 = 1;
                end
                if findedge(tempG, string(Last_dest), string(edge2_source)) > 0
                    isedge3 = 1;
                end
                if findedge(tempG, string(Last_dest), string(edge2_target)) > 0
                    isedge4 = 1;
                end

                if (isedge1 == 1) || (isedge2 == 1)|| (isedge3 == 1) || (isedge4 == 1)
                    over_ITER = over_ITER + 1;
                    continue
                end
                ttempG = tempG;
                ttempG = addedge(ttempG, string(First_dest), string(edge2_source));
                ttempG = addedge(ttempG, string(edge2_target), string(Last_dest));
                ttempG = rmedge(ttempG, string(edge2_source), string(edge2_target)); 
                [~, binsize] = conncomp(ttempG);
    
                if length(binsize) > 1
                    over_ITER = over_ITER + 1;
                    continue
                else
                    tempG = addedge(tempG, string(First_dest), string(edge2_source));
                    tempG = addedge(tempG, string(edge2_target), string(Last_dest));
                    tempG = rmedge(tempG, string(edge2_source), string(edge2_target));
                    changed_links = changed_links + 2;                    
                end
            else
                if over_ITER > 500
                    break
                end
                if Module_number(edge1_target) == Module_number(edge2_source)
                    continue
                end
                isedge2 = 0;
                isedge4 = 0;
                
                if findedge(tempG, string(edge1_target), string(edge2_target)) > 0
                    isedge2 = 1;
                end
                if findedge(tempG, string(edge1_target), string(edge2_source)) > 0
                    isedge4 = 1;
                end
                if (isedge2 == 1) && (isedge4 == 1)
                    continue
                end

                if isedge2 == 0          
                    ttempG = tempG;
                    ttempG = addedge(ttempG, string(edge1_target), string(edge2_target));
                    ttempG = rmedge(ttempG, string(edge2_source), string(edge2_target));              
                    [~, binsize] = conncomp(ttempG);

                    if length(binsize) > 1
                        over_ITER = over_ITER + 1;
                        continue
                    else
                        if (W - changed_links) == 3
                            Last_dest = edge2_source;
                        end
                        tempG = addedge(tempG, string(edge1_target), string(edge2_target));
                        tempG = rmedge(tempG, string(edge2_source), string(edge2_target));
                        edge1_target = edge2_source;
                        changed_links = changed_links + 1;                    
                    end
                    
                elseif isedge4 == 0
                    ttempG = tempG;
                    ttempG = addedge(ttempG, string(edge1_target), string(edge2_source));
                    ttempG = rmedge(ttempG, string(edge2_source), string(edge2_target));         
                    [~, binsize] = conncomp(ttempG);

                    if length(binsize) > 1
                        over_ITER = over_ITER + 1;
                        continue
                    else
                        if (W - changed_links) == 3
                            Last_dest = edge2_target;
                        end
                        tempG = addedge(tempG, string(edge1_target), string(edge2_source));
                        tempG = rmedge(tempG, string(edge2_source), string(edge2_target));
                        edge1_target = edge2_target;        
                        changed_links = changed_links + 1;                    
                    end
                end

                temp_intra(selected_idx) = [];
            end
        end
        
        
        tempA = adjacency(tempG);
        Q(j) = Q(j) + modularity(tempA, Module_number);
        %R(j) = R(j) + Robustness_BC(tempG);
        
        %R(j) = R(j) + Robustness_RF(tempG);

        Removelist = M_based_attack(tempG, Module_number);
        R(j) = R(j) + Robustness_MBA(tempG, Removelist);
        
        

    end
end

Q = Q./ITER;
R = R./ITER;
