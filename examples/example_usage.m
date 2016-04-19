simple_graphs;

input('Press any key to start graph draw.');
sgraphdraw(W1);

input('Press any key to start graph draw 2.');
sgraphdraw(A3bn);

input('Press any key to start main cluster demo.');
create_probablistic_cluster;
%sncut(A,2);
D=sum(A);
[IDX, XXc] = sncut(A(D>0,D>0),2);
fprintf('example cluster index.\n');
IDX
fprintf('example cluster vector XXc\n');
XXc

input('Press any key to continue cluster demo ... soft cluster.');
create_probablistic_cluster;
%sncut(A,2);
D=sum(A);
[IDX, XXc] = sncut(A(D>0,D>0),2, 'type','soft');
fprintf('example cluster index.\n');
IDX
fprintf('example cluster vector XXc\n');
XXc

input('Press any key to start graph draw with word labels.');
% example graph with words
W = [ ...            
    0    0.0442    0.0189    0.0150; ...
    0.0442         0    0.1585    0.0182; ...
    0.0189    0.1585         0    0.0097; ...
    0.0150    0.0182    0.0097         0; ...
    ];
example_words = { 'happy', 'abandon', 'abandonment', 'aggression' };
sgraphdraw(W, example_words);

input('Press any key to start interactive demo');
% example with prompt interface
ncutK_v5(A(D>0,D>0),5,5,0,0,1)