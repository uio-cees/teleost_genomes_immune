
import java.util.*;

/**
 * A node in the graph.
 *
 * @author Shalin Shah
 * Email: shah.shalin@gmail.com
 */
public class Node implements Comparable
{    
    public int value;
    public int degree = 0;
    public LinkedHashSet list = null;;
    public Node next = null;
    public Node previous = null;
    public int color = Constants.UNCOLORED;
    public int degreeSat = 0;
    public List possibleColors;
    public int colorCount = 0;
    
    /** Creates a new instance of Node */
    public Node(int v) 
    {
        value = v;
        list = new LinkedHashSet();
        possibleColors = new ArrayList();
    }
    
    /* 
     * compareTo method of the interface Comparable 
     */
    public int compareTo(Object o)
    {
        Node obj = (Node)o;
        
        if(this.value == obj.value)
        {
            return 0;
        }
        
        if(this.degree < obj.degree)
            return 1;
        else if(this.degree > obj.degree)
            return -1;
        else
            return 0;
    }
    
    /* 
     * Is the passed in object equal to this object? 
     */
    public boolean equals(Object obj)
    {
        Node o = (Node)obj;
        if(o.value == this.value)
            return true;
        
        return false;
    }
    
    /* Add a sibling */
    public void addEdge(int ev)
    {
        list.add(new Integer(ev));
    }
    
    /* 
     * Canonical String representation of this Node
     */
    public String toString()
    {
        return "Value=" + value + " Degree=" + degree;
    }
    
    public void colorNode(int col)
    {
        color = col;
    }
    
    public boolean isValidColor(Graph graph, int color)
    {
        Iterator it = list.iterator();
        while(it.hasNext())
        {
            int vertex = ((Integer)it.next()).intValue();
            Node node = graph.nodes[vertex];
            if(node.color == color)
            {   
                return false;
            }
        }
        
        return true;
    }
    
    public Node next()
    {
        return next;
    }
    
    public Node previous()
    {
        return previous;
    }
    
    public void computeDegreeSat(Graph graph)
    {
        degreeSat = 0;
        Iterator it = list.iterator();
        while(it.hasNext())
        {
            int vertex = ((Integer)it.next()).intValue();
            Node node = graph.nodes[vertex];
            if(node.color != Constants.UNCOLORED)
            {
                degreeSat++;
            }
        }
    }
 
    public int nextColor()
    {
        if(colorCount == possibleColors.size())
        {
            resetColorCount();
            color = Constants.UNCOLORED;
            return Constants.UNCOLORED;
        }
        
        int col = ((Integer)possibleColors.get(colorCount)).intValue();
        colorCount++;
        return col;
    }
    
    public void resetColorCount()
    {
        colorCount = 0;
        color = Constants.UNCOLORED;
    }
    
    public void computePossibleColors(Graph graph, int k)
    {
        possibleColors = new ArrayList();
        for(int i=1; i<=k; i++)
        {
            Iterator it = list.iterator();
            boolean flag = true;
            while(it.hasNext())
            {
                int vertex = ((Integer)it.next()).intValue();
                Node node = graph.nodes[vertex];
                if(node.color == i)
                {
                    flag = false;
                    break;
                }
            }
            
            if(flag)
            {
                possibleColors.add(new Integer(i));
            }
        }
    }
    
    public List findConflictingNodes(Graph graph)
    {
        List conflicts = new ArrayList();
        Iterator it = list.iterator();
        while(it.hasNext())
        {
            int vertex = ((Integer)it.next()).intValue();
            Node node = graph.nodes[vertex];
            if(node.color == this.color)
            {
                conflicts.add(node);
            }
        }
        
        return conflicts;
    }
    
    
    public int hashCode()
    {
        StringBuffer buffer = new StringBuffer("");
        buffer.append(value);
        buffer.append(" ");
        buffer.append(degree);
        return new String(buffer.toString()).hashCode();
    }
}