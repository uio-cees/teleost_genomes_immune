import java.util.*;

/**
 * Implementation of the DSatur heuristic for graph coloring.
 *
 * @author Shalin Shah
 * Email: shah.shalin@gmail.com
 */
public class GraphColoring {
    /** Creates a new instance of GraphColoring */
    public GraphColoring() {
    }
    
    public static void main(String [] args) throws Exception
    {
        System.err.println("Reading Graph...");
        /* Read the graph from Constants.FILE */
        Graph graph = GraphReader.readGraph();
        
        ArrayList<ArrayList<Node>> independentSets = new ArrayList<ArrayList<Node>>();
        
        /* Compute Clique */
        LinkedHashSet clique = MaxClique.computeClique(graph);
        
        /* Color the vertices of the clique with different colors */
        int k = clique.size();
        Iterator it = clique.iterator();
        int col = 1;        
        while(it.hasNext())
        {
            int vertex = ((Integer)it.next()).intValue();
            Node node = graph.nodes[vertex];
            node.colorNode(col);
            ArrayList<Node> ins = new ArrayList<Node>();
            ins.add(node);
            independentSets.add(ins);
            col++;
        }
            
        Node [] nodes = graph.nodes;
        PossibleColorsComparator comparator = new PossibleColorsComparator();
        TreeSet uncoloredNodes = new TreeSet(comparator);
        List listNodes = new ArrayList();
        LinkedHashSet flags = new LinkedHashSet();
        for(int i=0; i<nodes.length; i++)
        {
            Node node = nodes[i];
            if(node.color == Constants.UNCOLORED)
            {
                node.computeDegreeSat(graph);
                uncoloredNodes.add(node);
                listNodes.add(node);
                flags.add(new Integer(node.value));
            }
        }
        
        while(!uncoloredNodes.isEmpty())
        {
            Node node = (Node)uncoloredNodes.first();
            uncoloredNodes.remove(node);
            flags.remove(new Integer(node.value));
            //System.err.println(uncoloredNodes.size());
            
            final Integer[] idx = new Integer [k];
            final Integer[] siz = new Integer [k];
            for (int i = 0; i < k; i++) { 
            	idx[i] = i + 1;
            	siz[i] = independentSets.get(i).size();
            }
            //System.err.println(Arrays.toString(siz));
            Arrays.sort(idx, new Comparator<Integer>() {
                @Override public int compare(final Integer o1, final Integer o2) {
                    return siz[o1-1].compareTo(siz[o2-1]);
                }
            });
            
            for(int j=0; ;j++)
            {
            	int i = -1;
            	if (j < k ){
            		i = idx[j];
            	} else {
            		i = k + 1;
            		independentSets.add(new ArrayList<Node>());
            	}
            	//System.err.println(i);
                if(node.isValidColor(graph, i))
                {
                    node.colorNode(i);
                    independentSets.get(i-1).add(node);
                    LinkedHashSet list = node.list;
                    it = list.iterator();
                    while(it.hasNext())
                    {
                        int vertex = ((Integer)it.next()).intValue();
                        Node n = graph.nodes[vertex];
                        
                        if(uncoloredNodes.contains(n))
                        {
                            uncoloredNodes.remove(n);       
                            n.computeDegreeSat(graph);
                            uncoloredNodes.add(n);
                        }
                    }
                    
                    if(i > k)
                    {
                        k = i;
                    }
                    break;
                }
            }
        }
        
        System.err.println(k + " coloring found using DSatur.");
        //System.err.println("Applying Iterated Greedy Improvement...");
        
        //int [] colors = IteratedGreedy.iteratedGreedy(k, graph);
        int [] colors = new int[graph.nodes.length];
        
        for(int i=0; i<graph.nodes.length; i++)
        {
            colors[i] = graph.nodes[i].color;
        }
        
        int maxColor = -1;
        for(int i=0; i<colors.length; i++)
        {
            graph.nodes[i].color = colors[i];
            if(maxColor == -1)
            {
                maxColor = colors[i];
            }
            else if(maxColor < colors[i])
            {
                maxColor = colors[i];
            }
        }
        
/*        System.err.println("Applying Local Search...");
        colors = LocalSearch.localSearch(graph, maxColor);
        maxColor = -1;
        for(int i=0; i<colors.length; i++)
        {
            if(maxColor == -1)
            {        System.err.println("Applying Local Search...");
        colors = LocalSearch.localSearch(graph, maxColor);
        maxColor = -1;
        for(int i=0; i<colors.length; i++)
        {
            if(maxColor == -1)
            {
                maxColor = colors[i];
            }
            else if(maxColor < colors[i])
            {
                maxColor = colors[i];
            }
        }
                maxColor = colors[i];
            }
            else if(maxColor < colors[i])
            {
                maxColor = colors[i];
            }
        }*/
        
        System.err.println("Final Coloring of graph possible with " + maxColor + " colors.");
        System.err.println("Colors of Vertices: ");
        ArrayList<ArrayList<Integer>> res = new ArrayList<ArrayList<Integer>>();
        for(int i=0; i<maxColor; i++)                                                                                                            
        {                                                                                                                                        
            res.add(new ArrayList<Integer>());                                                                                                   
        }                                                                                                                                        
        for(int i=0; i<colors.length; i++)                                                                                                       
        {                                                                                                                                        
            System.err.print(colors[i] + " ");                                                                                                   
            res.get(colors[i]-1).add(i);                                                                                                         
        }                                                                                                                                        
        String ress = res.toString();                                                                                                            
        ress = ress.replace("]","|").replace("[","|").replace(",","");
        System.err.println(); 
        System.out.println(ress);   
    }
    
    public static class PossibleColorsComparator implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            Node node1 = (Node)o1;
            Node node2 = (Node)o2;
            
            if(node1.value == node2.value)
            {
                return 0;
            }
            
            if(node1.degreeSat < node2.degreeSat)
            {
                return 1;
            }
            else
            {
                return -1;
            }
        }
    }
}

