
import java.util.*;

/**
 * Backtracking solver for the graph coloring problem.
 * Spends a specific amount of time (Constants.TIME) on trying to color the graph
 * using k colors after which it increments k by 1.
 * If k is known, it can be populated (Constants.KNOWN_K).
 * The time it spends on each k before incrementing it is very important for the performance of the
 * algorithm. For large graphs, the time should be increased to at least 10000 ms.
 * For small dense graphs, 2000 ms should work.
 *
 * @author Shalin Shah
 * Email: shah.shalin@gmail.com
 */
public class Backtracking {
    
    /** Creates a new instance of GraphColoring */
    public Backtracking() {
    }
    
    public static void main(String [] args) throws Exception
    {
        System.err.println("Reading Graph...");
        /* Read the graph from Constants.FILE */
        Graph graph = GraphReader.readGraph();
        
        /* Compute Clique */
        LinkedHashSet clique = MaxClique.computeClique(graph);
        
        /* Color the vertices of the clique with different colors */
        int k = clique.size();
        if(Constants.KNOWN_K > 0)
            k = Constants.KNOWN_K;
        
        Iterator it = clique.iterator();
        int col = 0;
        while(it.hasNext())
        {
            int vertex = ((Integer)it.next()).intValue();
            Node node = graph.nodes[vertex];
            node.colorNode(col);
            col++;
        }
        
        /* Try to color the nodes with increasing values of k */
        
        while(true)
        {
            List uncoloredNodes = new ArrayList();
            Node [] nodes = graph.nodes;
            
            for(int i=0; i<nodes.length; i++)
            {
                Node node = nodes[i];
                if(node.color == Constants.UNCOLORED)
                {
                    node.computePossibleColors(graph, k);
                    uncoloredNodes.add(node);
                }
            }
            
            PossibleColorsComparator comparator = new PossibleColorsComparator();
            Collections.sort(uncoloredNodes, comparator);
            Node previous = null;
            for(int i=0; i<uncoloredNodes.size();i++)
            {
                Node node = (Node)uncoloredNodes.get(i);
                if(previous == null)
                {
                    previous = node;
                    node.previous = null;
                }
                else
                {
                    previous.next = node;
                    node.previous = previous;
                    previous = node;
                }
            }
            Node head = (Node)uncoloredNodes.get(0);
            
            /* Try to color the uncolored nodes using k colors */
            Node current = (Node)uncoloredNodes.get(0);
            Node last = (Node)uncoloredNodes.get(uncoloredNodes.size()-1);
            boolean solved = false;
            System.err.println("Trying to color graph with " + k + " colors...");
            long startTime = System.currentTimeMillis();
            while(true)
            {
                long endTime = System.currentTimeMillis();
                if((endTime-startTime) > Constants.TIME)
                {
                    solved = false;
                    break;
                }
                
                if(current == last)
                {
                    int nextColor = current.nextColor();
                    while(nextColor != Constants.UNCOLORED && !current.isValidColor(graph, nextColor))
                    {
                        nextColor = current.nextColor();
                    }
                    
                    if(nextColor == Constants.UNCOLORED)
                    {
                        current.resetColorCount();
                        current = current.previous;
                    }    
                    else
                    {
                        /* Solved */
                        current.colorNode(nextColor);
                        System.err.println("\n" + k + " coloring found! Exiting...");
                        solved = true;
                        break;
                    }   
                }
                else
                {
                    int nextColor = current.nextColor();
                    while(nextColor != Constants.UNCOLORED && !current.isValidColor(graph, nextColor))
                    {
                        nextColor = current.nextColor();
                    }
                    
                    if(nextColor == Constants.UNCOLORED)
                    {
                        current.resetColorCount();
                        current = current.previous;
                        if(current == null)
                        {
                            solved = false;
                            break;
                        }
                    }
                    else
                    {
                        current.colorNode(nextColor);
                        current = current.next();
                    }
                }
            }
            
            if(solved)
            {
                break;
            }
            else
            {
                k++;
            }
        }
    }
    
    public static class PossibleColorsComparator implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            Node node1 = (Node)o1;
            Node node2 = (Node)o2;
            if(node1.possibleColors.size() < node2.possibleColors.size())
            {
                return 1;
            }
            else if(node1.possibleColors.size() > node2.possibleColors.size())
            {
                return -1;
            }
            else
            {
                return 0;
            }
        }
    }
}
