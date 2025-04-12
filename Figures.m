function main()
    [C,X,I] = loadDS();
    sigma = 1;
    [Sc, Ss, Si] = sim_matrices(C, X, I, sigma);
    plot_scatt(C, "Circle dataset")
    plot_scatt(X, "Spiral dataset")
    plot_scatt(I, "Iris dataset")

    [Wc, Dc, Ws10, Ds10, Ws40, Ds40, Wi, Di] = construct_graphs(C, Sc, X, Ss, I, Si);

    plotLaplacians(Wc, Dc, Ws10, Ds10, Ws40, Ds40, Wi, Di);   
    
    end

function [C,X, I] = loadDS()
    datasetStructC = load("Circle.mat");
    C = datasetStructC.X;

    datasetStructS = load("Spiral.mat");
    X = datasetStructS.X;
    X = X(:, 1:2);  
    
    datasetStructI = load("Iris.mat");
    I = datasetStructI.X;
    I= I(:, 1:3);  
end

function [Sc, Ss, Si] = sim_matrices(C, X, I, sigma)
    Dic = squareform(pdist(C));
    Dis = squareform(pdist(X));
    Dii = squareform(pdist(I));
    
    Dic = sparse(Dic);
    Dis = sparse(Dis);
    Dii= sparse(Dii);

    Sc = exp(-(Dic.^2 / (2 * sigma^2)));
    Ss = exp(-(Dis.^2 / (2 * sigma^2)));
    Si = exp(-(Dii.^2 / (2 * sigma^2)));
    
    Sc(eye(size(Sc)) == 1) = 0;
    Ss(eye(size(Ss)) == 1) = 0;
    Si(eye(size(Si)) == 1) = 0;
    
    Sc(Sc < 10e-7) = 0;
    Ss(Ss < 10e-7) = 0;
    Si(Si < 10e-7) = 0;


end
function [Wc,Dc,Ws10,Ds10,Ws40,Ds40,Wi,Di] = construct_graphs(C, Sc, X, Ss, I, Si)
    [Wc, Dc] = construct_graph_matrix(C, Sc,10);
    [Ws10, Ds10] = construct_graph_matrix(X, Ss,10);
    [Ws40, Ds40] = construct_graph_matrix(X, Ss,40);
    [Wi, Di] = construct_graph_matrix(I, Si,10); 
       
    Gc = graph(Wc);
    figure;
    plot(Gc, "Xdata", C(:, 1), "Ydata",C(:, 2)); 
    title('10-Nearest Neighbors Similarity Graph');
    subtitle(sprintf('Number of connected components: %d', max(conncomp(Gc))), 'Units', 'normalized');
  end

function [W, D] = construct_graph_matrix(data, similarity,k)
    nPoints = size(data, 1);
    W = zeros(nPoints, nPoints);
    
    for i = 1:nPoints
        [~, indices] = sort(similarity(i, :), 'descend');
        neighbors = indices(1:k);  
        W(i, neighbors) = similarity(i, neighbors);
    end
    
    W = max(W, W'); 
    W=sparse(W);
    D = diag(sum(W, 2));
    D=sparse(D);
end

function plot_scatt(data, titleText)
    figure;

    if strcmp(titleText, "Iris dataset")
        view(3); 
        scatter3(data(:, 1), data(:, 2), data(:, 3), 16, 'o', 'filled'); 
        xlabel('X1');
        ylabel('X2');
        zlabel('X3');
        title('Iris dataset');
    else
        scatter(data(:, 1), data(:, 2), 16, 'o', 'filled');
        title(sprintf('Plot of %s', titleText));
        xlabel('X1');
        ylabel('X2');
    end
end
function plotLaplacians(Wc, Dc, Ws10, Ds10, Ws40, Ds40, Wi, Di)
    figure;
    Lc=Dc-Wc;
    spy(Lc);
    title(sprintf('Sparsity pattern graph of L'));
    subtitle(sprintf('For k=10, Circle dataset'));

    figure;
    Li=Di-Wi;
    spy(Li);
    title(sprintf('Sparsity pattern graph of L'));
    subtitle(sprintf('For k=10, Iris dataset'));


    figure;
    Ls10=Ds10-Ws10;
    spy(Ls10);
    title(sprintf('Sparsity pattern graph of L'));
    subtitle(sprintf('For k=10, Spiral dataset'));

    
    figure;
    Ls40=Ds40-Ws40;
    spy(Ls40);
    title(sprintf('Sparsity pattern graph of L'));
    subtitle(sprintf('For k=40, Spiral dataset'));



end
