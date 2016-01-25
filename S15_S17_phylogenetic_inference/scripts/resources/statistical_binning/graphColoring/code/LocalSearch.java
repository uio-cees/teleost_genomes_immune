
import java.util.*;

/**
 * Min-Conflicts local search heuristic
 * @author Shalin Shah
 * Email: shah.shalin@gmail.com
 */
public class LocalSearch 
{    
    /* 
     * Take a valid coloring and try to improve the coloring
     * by reducing one color and then trying to remove the conflicts
     * using the min-conflicts heuristic
     */
    public static int [] localSearch(Graph graph, int k)
    {
        long startTime = System.currentTimeMillis();
        int [] coloring = new int[graph.nodes.length];
        for(int i=0; i<coloring.length;  i++)
        {
            coloring[i] = graph.nodes[i].color;
        }
        
        for(int p=0; p<Constants.LOCAL_SEARCH_ITERATIONS;p++)
        {
            long endTime = System.currentTimeMillis();
            if(endTime - startTime > Constants.LOCAL_SEARCH_MAX_TIME)
            {
                break;
            }
            
            int maxColor = k-1;
            for(int i=0; i<graph.nodes.length; i++)
            {
                Node node = graph.nodes[i];
                if(node.color == k)
                {
                    int c = ((int)(Math.random()*(double)maxColor)) + 1;
                    node.color = c;
                }   
            }
            
            boolean flag = false;
            for(int n=0; n<1000; n++)
            {
                List conflicts = new ArrayList();
                for(int i=0; i<graph.nodes.length; i++)
                {
                    Node node = graph.nodes[i];
                    if(!node.isValidColor(graph, node.color))
                    {
                        conflicts.add(node);
                    }
                }
                //System.err.println(conflicts.size());
                
                if(conflicts.size() == 0)
                {
                    System.err.println("Found Better Coloring - " + maxColor);
                    flag = true;
                    for(int i=0; i<coloring.length; i++)
                    {
                        coloring[i] = graph.nodes[i].color;
                    }
                    k--;
                    break;
                }
                else
                {
                    changeColorsRandomly(conflicts, maxColor);
                    conflicts = new ArrayList();
                    for(int i=0; i<graph.nodes.length;i++)
                    {
                        Node node = graph.nodes[i];
                        if(!node.isValidColor(graph, node.color))
                        {
                            conflicts.add(node);
                        }
                    }
                }
                
                for(int i=0; i<1000; i++) 
                {
                    if(conflicts.size() == 0)
                    {
                        //System.err.println("Better Coloring Found!");
                        break;
                    }
                    
                    int rand = (int)(Math.random() * conflicts.size());
                    Node node = (Node)conflicts.get(rand);
                    int bestcolor = -1;
                    int bestconflicts = -1;
                    for(int c=1; c<=maxColor; c++) 
                    {
                        node.color = c;
                        List con = node.findConflictingNodes(graph);
                        
                        if(con.size() == 0) 
                        {
                            bestcolor = c;
                            conflicts.remove(node);
                            break;
                        }
                        
                        if(bestcolor == -1) 
                        {
                            bestcolor = c;
                            bestconflicts = con.size();
                        } 
                        else 
                        {
                            if(bestconflicts > con.size()) 
                            {
                                bestconflicts = con.size();
                                bestcolor = c;
                            }
                        }
                    }
                    
                    node.color = bestcolor;
                }
            }
            
            if(!flag)
            {
                break;
            }
        }
        
        return coloring;
    }
    
    public static void changeColorsRandomly(List conflicts, int k)
    {
        Iterator it = conflicts.iterator();
        while(it.hasNext())
        {
            Node node = (Node)it.next();
            int color = ((int)(Math.random()*k)) + 1;
            node.color = color;
        }
    }
}
