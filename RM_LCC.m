
W_list = [0.987362, 0.937368, 0.874875, 0.624903, 0.374931];
N = numnodes(G);
ITER = 100;
LCC_IB = zeros(length(W_list),N);
SLCC_IB = zeros(length(W_list),N);

LCC_RF = zeros(length(W_list),N);
SLCC_RF = zeros(length(W_list),N);
A = adjacency(G);

%%
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
        temp_inter = Inter_edges;
        W = round(W_list(j)*length(Inter_edges));
        changed_links = 0;
        First_selected_idx = randperm(length(temp_inter),1);
        First_selected_edges = temp_inter(First_selected_idx);
        edge1 = Edgelist{First_selected_edges,['EndNodes']};
        edge1_source = str2double(edge1(1));
        edge1_target = str2double(edge1(2));
        temp_inter(First_selected_idx) = [];
        tempG = rmedge(tempG, string(edge1_source), string(edge1_target));
        First_dest = edge1_source;

        over_ITER = 0;
        isnoedge = 0;

        while changed_links < W
            isskip = 0;
            selected_idx = randperm(length(temp_inter),1);
            selected_edges = temp_inter(selected_idx);
            edge2 = Edgelist{selected_edges,['EndNodes']};
            edge2_source = str2double(edge2(1));
            edge2_target = str2double(edge2(2));

            if (W - changed_links) < 3
                if over_ITER > 500
                    break
                end
                if Module_number(First_dest) ~= Module_number(edge2_source)
                    over_ITER = over_ITER + 1;
                    continue
                end
                if Module_number(edge2_target) ~= Module_number(Last_dest)
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
                isModSame2 = 0;
                isModSame4 = 0;
                
                if Module_number(edge1_target) == Module_number(edge2_source)
                    isModSame2 = 1;
                end
                if Module_number(edge1_target) == Module_number(edge2_target)
                    isModSame4 = 1;
                end
                
                isedge2 = 0;
                isedge4 = 0;
                
                if findedge(tempG, string(edge1_target), string(edge2_source)) > 0
                    isedge2 = 1;
                end
                
                if findedge(tempG, string(edge1_target), string(edge2_target)) > 0
                    isedge4 = 1;
                end
                
                isOK2 = 0;
                isOK4 = 0;
                
                if (isModSame2 == 1) && (isedge2 == 0)
                    isOK2 = 1;
                end
                if (isModSame4 == 1) && (isedge4 == 0)
                    isOK4 = 1;
                end
                
                if (isOK2 == 0) && (isOK4 == 0)
                    isnoedge = isnoedge + 1;
                    if isnoedge > (length(temp_inter)*2)
                        if isedge2 == 0
                            isOK2 = 1;
                            isskip = 1;
                        end
                        if isedge4 == 0
                            isOK4 = 1;
                            isskip = 1;
                        end
                    else
                        continue
                    end
                end


                if isOK2 == 1          
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
                        if isskip == 0
                            changed_links = changed_links + 1;
                        end
                    end
                    
                elseif (isOK2 == 0) && (isOK4 == 1)
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
                        if isskip == 0
                            changed_links = changed_links + 1;
                        end
                    end
                end

                temp_inter(selected_idx) = [];
                isnoedge = 0;
            end
        end
        
        [~, Removelist] = HBCA_Non_renewal(tempG, 1);
        ttempG = tempG;
        for i = 1:(length(Removelist)-1)
            ttempG = rmnode(ttempG, string(Removelist(i)));
            [~, binsize] = conncomp(ttempG);
            FandS = maxk(binsize,2);
            if length(FandS) > 1
                LCC_IB(j,i) = LCC_IB(j,i) + FandS(1)/N;
                SLCC_IB(j,i) = SLCC_IB(j,i) + FandS(2)/N;
            else
                LCC_IB(j,i) = LCC_IB(j,i) + FandS(1)/N;
                SLCC_IB(j,i) = SLCC_IB(j,i) + 0;
            end
        end
        
        Removelist = randperm(N);
        ttempG = tempG;
        for i = 1:(length(Removelist)-1)
            ttempG = rmnode(ttempG, string(Removelist(i)));
            [~, binsize] = conncomp(ttempG);
            FandS = maxk(binsize,2);
            if length(FandS) > 1
                LCC_RF(j,i) = LCC_RF(j,i) + FandS(1)/N;
                SLCC_RF(j,i) = SLCC_RF(j,i) + FandS(2)/N;
            else
                LCC_RF(j,i) = LCC_RF(j,i) + FandS(1)/N;
                SLCC_RF(j,i) = SLCC_RF(j,i) + 0;
            end
        end
        
    end
end
LCC_IB = LCC_IB./ITER;
SLCC_IB = SLCC_IB./ITER;
LCC_RF = LCC_RF./ITER;
SLCC_RF = SLCC_RF./ITER;
