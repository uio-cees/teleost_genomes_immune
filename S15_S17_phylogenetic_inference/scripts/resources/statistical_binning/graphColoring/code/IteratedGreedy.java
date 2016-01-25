
import java.util.*;

/**
 * Iterated Greedy improvement procedure.
 *
 * @author Shalin Shah
 * Email: shah.shalin@gmail.com
 */
public class IteratedGreedy {
    
    /** Creates a new instance of IteratedGreedy */
    public IteratedGreedy() {
    }
    
    public static int [] iteratedGreedy(int k, Graph graph)
    {
        int [] colors = new int[graph.nodes.length];
        
        for(int i=0; i<graph.nodes.length; i++)
        {
            colors[i] = graph.nodes[i].color;
        }
        
        for(int i=0; i<Constants.ITERATIONS; i++)
        {
            double rand = Math.random();
            Node [] order = null;
            double left = 0;
            double right = 0;
            right+=Constants.RANDOM;
            boolean chosen = false;
            if(rand < right && rand > left && !chosen)
            {
                order = randomOrdering(graph.nodes);
                chosen = true;
            }
            
            left = right;
            right+=Constants.DECREASING;
            if(rand < right && rand > left && !chosen)
            {
                order = decreasingOrdering(graph.nodes);
                chosen = true;
            }
            
            left = right;
            right+=Constants.INCREASING;
            
            if(rand < right && rand > left && !chosen)
            {
                order = increasingOrdering(graph.nodes);
                chosen = true;
            }
            
            left = right;
            right+=Constants.REVERSE;
            
            if(rand < right && rand > left && !chosen)
            {
                order = reverseOrdering(graph.nodes);
                chosen = true;
            }
            
            int maxColor = -1;
            for(int n=0; n<order.length; n++)
            {
                Node node = order[n];
                node.color = Constants.UNCOLORED;
            }
            
            for(int n=0; n<order.length; n++)    
            {
                Node node = order[n];
                for(int c=1;;c++)
                {
                    if(node.isValidColor(graph, c))
                    {
                        node.colorNode(c);
                        if(maxColor == -1)
                        {
                            maxColor = c;
                        }
                        else if(maxColor < c)
                        {
                            maxColor = c;
                        }
                        break;
                    }
                }
            }
            
            if(maxColor < k)
            {
                System.err.println("Found Better Coloring - " + maxColor);
                k = maxColor;
                for(int j=0; j<colors.length; j++)
                {
                    colors[j] = graph.nodes[j].color;
                }
            }
        }
        
        return colors;
    }
    
    public static Node [] randomOrdering(Node [] order)
    {
        List ns = new ArrayList(Arrays.asList(order));
        Collections.shuffle(ns);
        Node [] nodes = new Node[order.length];
        for(int i=0; i<order.length; i++)
        {
            nodes[i] = (Node)ns.get(i);
        }
        
        return nodes;
    }
    
    public static Node [] decreasingOrdering(Node [] order)
    {
        int k = -1;
        for(int i=0; i<order.length; i++)
        {
            if(k == -1)
            {
                k = order[i].color;
            }
            else if(k < order[i].color)
            {
                k = order[i].color;
            }
        }
        
        ColorClass [] classes = new ColorClass[k];
        for(int i=0; i<order.length; i++)
        {
            int color = order[i].color;
            if(classes[color-1] == null)
            {
                classes[color-1] = new ColorClass();
            }
            classes[color-1].nodes.add(order[i]);
        }
        
        Arrays.sort(classes, new DecreasingComparator());
        for(int i=0; i<classes.length; i++)
        {
            Collections.reverse(classes[i].nodes);
        }
        
        Node [] nodes = new Node[order.length];
        int count = 0;
        for(int i=0; i<classes.length; i++)
        {
            ColorClass cls = classes[i];
            List nds = cls.nodes;
            for(Iterator it = nds.iterator(); it.hasNext();)
            {
                nodes[count] = (Node)it.next();
                count++;
            }
        }
        
        return nodes;
    }
    
    public static Node [] increasingOrdering(Node [] order)
    {
        int k = -1;
        for(int i=0; i<order.length; i++)
        {
            if(k == -1)
            {
                k = order[i].color;
            }
            else if(k < order[i].color)
            {
                k = order[i].color;
            }
        }
        
        ColorClass [] classes = new ColorClass[k];
        for(int i=0; i<order.length; i++)
        {
            int color = order[i].color;
            if(classes[color-1] == null)
            {
                classes[color-1] = new ColorClass();
            }
            classes[color-1].nodes.add(order[i]);
        }
        
        Arrays.sort(classes, new IncreasingComparator());
        for(int i=0; i<classes.length; i++)
        {
            Collections.reverse(classes[i].nodes);
        }
        
        Node [] nodes = new Node[order.length];
        int count = 0;
        for(int i=0; i<classes.length; i++)
        {
            ColorClass cls = classes[i];
            List nds = cls.nodes;
            for(Iterator it = nds.iterator(); it.hasNext();)
            {
                nodes[count] = (Node)it.next();
                count++;
            }
        }
        
        return nodes;
    }
    
    public static Node [] reverseOrdering(Node [] order)
    {
        List ns = new ArrayList(Arrays.asList(order));
        Collections.reverse(ns);
        Node [] nodes = new Node[order.length];
        for(int i=0; i<nodes.length; i++)
        {
            nodes[i] = (Node)ns.get(i);
        }
        
        return nodes;
    }
    
    public static class ColorClass
    {
        public List nodes = new ArrayList();
        public int color = -1;
    }
    
    public static class DecreasingComparator implements Comparator
    {
        public int compare(Object obj1, Object obj2)
        {
            ColorClass c1 = (ColorClass)obj1;
            ColorClass c2 = (ColorClass)obj2;
            
            if(c1.nodes.size() > c2.nodes.size())
            {
                return -1;
            }
            else if(c1.nodes.size() < c2.nodes.size())
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }
    
    public static class IncreasingComparator implements Comparator
    {
        public int compare(Object obj1, Object obj2)
        {
            ColorClass c1 = (ColorClass)obj1;
            ColorClass c2 = (ColorClass)obj2;
            
            if(c1.nodes.size() > c2.nodes.size())
            {
                return 1;
            }
            else if(c1.nodes.size() < c2.nodes.size())
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
