
Heuristic for Graph Coloring
----------------------------

Implementation of the DSatur heuristic for graph coloring. The heuristic follows
the following steps:
	
	1) Compute a clique (maximum is good)
	2) Color the clique
	3) Sort the rest of the vertices in non-increasing order of the 
	degree of saturation
	4) Color the vertices in the order given by 3. Also, when a vertex is 
	colored, change the degree of saturation of the neighbouring vertices
	so that the order of coloring changes
	5) Improve the coloring using Iterated Greedy
	6) Improve the coloring using min-conflicts local search
	7) Report the coloring

Author: Shalin Shah
Email: shah.shalin@gmail.com