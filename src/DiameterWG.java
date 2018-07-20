import java.util.List;
import java.io.File;
import java.io.PrintStream;
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.Arrays;
import java.util.ArrayList;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.stat.StatUtils;
import java.util.Scanner;


public class DiameterWG {
  ImmutableGraph G; // WebGraph object
  int n; // number of nodes
  double[] score;
  Integer[] index;
  int diameter;


  public DiameterWG(String basename) throws Exception {
    G = ImmutableGraph.loadMapped(basename);
    n = G.numNodes();
    score= new double[n];
    Arrays.parallelSetAll(score, (i) -> 0);
    index = new Integer[n];
    Arrays.parallelSetAll(index, (i) -> i);
  }

 // compute the dimeter
  public int scoreCompute() {

    score = Arrays.stream(index)
                  .parallel()
                  .mapToDouble( u -> {
                  int source = u;
                  System.out.println("Start visiting node " + source);
                  ImmutableGraph H = G.copy();
                  // queue stores the visited nodes in BFS order
                  List<Integer> queue = new ArrayList<Integer>();
                  // breakPoint stores start positions (index in queue) of each BFS layer
                  // For example:
                  // breakPoint[0]=0  -->     *         depth = 0
                  // breakPoint[1]=1  -->    * *        depth = 1
                  // breakPoint[2]=3  -->   * * *       depth = 2
                  // breakPoint[3]=6  -->  * * * *      depth = 3
                  // breakPoint[4]=10 --> * * * * *     depth = 4
                  //
                  // breakPoint[0] = 0 (index (in queue) of the root node),
                  // breakPoint[1] = 1 (index (in queue) of the first node of layer 1 (depth = 1))
                  // We can think this number as the number of nodes above this layer.
                  // breakPoint[2] = 3 (index (in queue) of the first node of layer 2 (depth = 2))
                  // or the number of nodes above this layer is 3.
                  // breakPoint[k] = (index (in queue) of the first node of layer k (depth = k))
                  // The last element of breakPoint will be the size of queue which is n.
                  List<Integer> breakPoint = new ArrayList<Integer>();
                  // distance[v] stores the distance from source to v
                  double[] distance = new double[n];
                  Arrays.parallelSetAll(distance, (i) -> -1);
                  distance[source] = 0;
                  // use BFS to find shortest paths
                  queue.clear();
                  queue.add(source);
                  breakPoint.clear();
                  breakPoint.add(0);

                  int depth; // depth of the BFS spanning tree
                  for(depth = 0; queue.size() != breakPoint.get(breakPoint.size() - 1); depth++) {
                    breakPoint.add(queue.size());
                    int start = breakPoint.get(depth);
                    int end = breakPoint.get(depth + 1);

                    for(int i = start; i < end; i++) {
                      int v = queue.get(i);
                      for(int w : H.successorArray(v)) {
                        if (distance[w] == -1) {
                          distance[w] = depth + 1;
                          queue.add(w);
                        }
                      }
                    }
                  }
                 System.out.println("Ending BFS in node " + source);
                 // finding the largest shortest path
                  return StatUtils.max(distance);})
                  .toArray();
    // finding the diameter
    diameter = (int) StatUtils.max(score);
    return diameter;
  }


  public static void main(String[] args) throws Exception {
    System.out.println("Please enter graph's basename.");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();

		long startTime = System.currentTimeMillis();

    // define extra path strings
    String inputPath = "data/";
    String outputPath = "/output/diameter/";

    PrintStream ps2 = new PrintStream(new File(inputPath+basename+"/"+outputPath+basename+"_diameter_WG_timing.txt"));

		System.out.println("Starting " + basename);
    DiameterWG d = new DiameterWG(inputPath+basename+"/"+basename);

    System.out.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);
    ps2.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);

    long t2 = System.currentTimeMillis();
	  int res = d.scoreCompute();
    System.out.println("Diameter = " + res);
    ps2.println("Diameter = " + res);

    System.out.println("Time elapsed (sec) for diameter computing = " + (System.currentTimeMillis() - t2)/1000.0);
    ps2.println("Time elapsed (sec) for score computing = " + (System.currentTimeMillis() - t2)/1000.0);


    System.out.println("Total time elapsed (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);
    ps2.println("Total time elapsed (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);

	}
}
