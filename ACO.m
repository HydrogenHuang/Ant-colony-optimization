function ACO
%Figure out TSP problem by Ant colony optimization

m=3;
n=4;
alpha=1;
beta=2;
rho=0.5;
NC_max = 100;%the maximun of iteration
oo=eps;%no distance
D = [oo 3 1 2;
     3 oo 5 4;
     1 5 oo 2;
     2 4 2 oo];%the distance between two cities.

%Initailize
Eta = 1./D;%heuristic factor
Tau = 0.3.*(ones(n,n)-eye(n,n));%pheromone matrix
R_best = zeros(NC_max,n);%the best route of each generation
L_best = inf.*ones(NC_max,1);%the length of the best route

NC = 1;%iterator
while NC<=NC_max
    Tabu = zeros(m,n);%record the trade
    %Random select places
    Tabu(:,1) = ceil(rand(m,1).*10./n);
    %Random select next city
    for j=2:n
        for i=1:m
            visited = Tabu(i,1:(j-1));%record the visited cities
            screen = zeros(1,(n-j+1));%cities to be screened
            
            %Find the unvisited cities
            idx = 1;
            for k=1:n
                if isempty(find(visited==k, 1))
                    screen(idx) = k;
                    idx = idx+1;
                end
            end
            
            %Calculate the possibility of each unvisited cities
            P = screen;
            for k=1:length(screen)
                P(k) = (Tau(visited(end),screen(k))^alpha)*(Eta(visited(end),screen(k))^beta);
            end
            P = P./sum(P);
            
            %Select the next city by probability
            Pcum = cumsum(P);
            select = find(Pcum>=rand);
            Tabu(i,j) = screen(select(1));
        end
    end
    
    %
    %
    
    %Record the best route
    L = zeros(m,1);%length
    for i=1:m
        R = Tabu(i,:);%route
        for j=1:(n-1)
            L(i) = L(i)+D(R(j),R(j+1));
        end
        L(i) = L(i)+D(R(n),R(1));
    end
    [L_best(NC),idx] = min(L);
    R_best(NC,:) = Tabu(idx,:);
    NC = NC+1;
    
    %Renew the heuristic factor
    Delta = zeros(n,n);
    for i=1:m
        for j=1:(n-1)
            Delta(Tabu(i,j),Tabu(i,j+1)) = Delta(Tabu(i,j),Tabu(i,j+1))+1/L(i);
        end
        Delta(Tabu(i,n),Tabu(i,1)) = Delta(Tabu(i,n),Tabu(i,1))+1/L(i);
    end
    Tau = (1-rho).*Tau+Delta;
end

%Print out the outcome
[Shortest_Length,idx] = min(L_best);
Shortest_Route = R_best(idx,:);%the best route
fprintf('The shortest route is: ');
for i=1:n
    fprintf('%d -> ',Shortest_Route(i));
end
fprintf('%d',Shortest_Route(1));
fprintf('\nThe length of this route is %d\n',Shortest_Length);
