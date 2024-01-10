function R = Robustness(G,loop)
n = numnodes(G);
%Sq_N = zeros(1,n);
R = zeros(1,loop);
for k = 1:loop
    S_q = zeros(1,n);
    temp = G;
    for i = 1:(n-1)
        Names = string(temp.Nodes.Name);
        C = zeros(2,length(Names));
        C(1,:) = Names';
        C(2,:) = centrality(temp, 'degree');
        d_max = max(C(2,:));
        d_max_index = find(C(2,:) == d_max);
        if size(d_max_index,2) > 1
            d_max_index = datasample(d_max_index,1);
        end
        temp = rmnode(temp,string(C(1,d_max_index)));
        [~, Cluster_size] = conncomp(temp);
        if isempty(Cluster_size) == 1
            Cluster_size = 0;
        end
        LCC = max(Cluster_size);
        S_q(i) = LCC/n;
    end
    %Sq_N = Sq_N+S_q/n;
    R(k) = sum(S_q)/(n-1);
end
%Sq_N = Sq_N/loop;
R = sum(R)/loop;

end

