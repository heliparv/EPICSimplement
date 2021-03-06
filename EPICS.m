%Number of species
%number of species from Gould et al (as referenced by Ansari et al)
n = 5;
%below number of species in study by Ansari et al
%n=8;

%abundance in monoculture
%below abundances from Gould et al (as referenced by Ansari et al)
abundance_in_monoc = [177827.941003893, 422668.614265604, 143301.257023696, 237137.370566166, 143301.257023696];
%below abundances from Ansari et al
%abundance_in_monoc = [2884030259.24548, 41686.9890738402, 1258925.78380719, 436515.384984744, 8709663.37151069, 199526.324505624, 346737.031168192, 2137962.73530952]

%abundances in leave-one-out cultures
%below abundances from Gould et al (as referenced by Ansari et al)
abundance_in_loo = [0.0, 396525.520961999, 123647.743095677, 170548.611166451, 115120.312537355; 254186.34555361, 0.0, 21405.1659413566, 69566.7893094091, 141809.224361488; 202495.347480754, 205986.64657525, 0.0, 104738.972834873, 136160.664685335; 242561.433240061, 272415.148100376, 59707.4297206303, 0.0,  123146.5737988; 255231.590042967,  241558.469147808,  195981.399497278,  173192.864672013, 0.0];
%below abundances from Ansari et al
%abundance_in_loo = [0, 91201.08394, 16595869.07, 831763.7711, 1174.897555, 36307.80548, 4897.788194, 11481.53621; 8317637.711, 0, 51286138.4, 9549925.86, 263026.7992, 16595.86907, 12022.64435, 2290.867653; 1047128.548, 48977.88194, 0, 35481338.92, 95499.2586, 46773.51413, 23988.32919, 331131.1215; 13.18256739, 70794.57844, 19952623.15, 0, 3162.27766, 18620.87137, 2630.267992, 758.577575; 281838293.1, 27542.28703, 1995262.315, 3235936.569, 0, 288.4031503, 4365.158322, 5623413.252; 467735.1413, 81283.05162, 186208713.7, 16595869.07, 95499.2586, 0, 15848.93192, 630.9573445; 20417.37945, 89125.09381, 407380.2778, 50118.72336, 22908.67653, 123.0268771, 0, 9549.92586; 218776.1624, 79432.82347, 218776162.4, 2630267.992, 12882.49552, 3019.95172, 15848.93192, 0];

%Creates matrix A for solving Ax=b, consists of species abundances in monoculture and l-o-o
abundance_matrix = zeros(n*n,n*n);
row = 1;
for i=1:n
    for j=1:n
        if i==j
            abundance_matrix(row,n*i-n+i) = abundance_in_monoc(i);
        else
            abundance_matrix(row,n*j-n+1:n*j) = abundance_in_loo(i,1:n);
        end
        row = row+1;
    end
end

%solves Ax=b creating a vector of effective pairwise interactions, b is a vector of -1 as growth rates are set to 1
interactions_vector = abundance_matrix\(-1*ones(n*n,1));

%creates matrix of effective pairwise interactions to be used as A in Ax=b
interactions_matrix = (reshape(interactions_vector,n,n))'

%normalizes interaction matrix by self-interactions like done in EPICS paper
normalized_interactions_matrix = interactions_matrix;
for i=1:n
    for j=1:n
        normalized_interactions_matrix(i,j) = normalized_interactions_matrix(i,j)/interactions_matrix(j,j);
    end
end
normalized_interactions_matrix

#calculates prediction of abundances in n-member community by solving Ax=b where A is effective pairwise interactions and b set to -1
abundance_n_member_community = (interactions_matrix\(-1*ones(n,1)))'