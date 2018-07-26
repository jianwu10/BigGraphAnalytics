/*
*  Effective diameter estimation using HyperLogLog counter
*/
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.ArrayList;
import java.util.Scanner;
import it.unimi.dsi.webgraph.algo.HyperBall;
import it.unimi.dsi.webgraph.algo.NeighbourhoodFunction;

public class DiameterHyperLogLog {
  public static void main(String[] args) throws Exception {
    System.out.println("Please enter graph's base name");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();

    System.out.println("Please enter the K value (e.g., 16, 32, 64,...)");
    int k_value = sc.nextInt();

    System.out.println("Please enter the max iteration number");
    int maxIter = sc.nextInt();

    long t1 = System.currentTimeMillis();

    ImmutableGraph G = ImmutableGraph.loadMapped("data/"+basename+"/"+basename);
    HyperBall hyperball = new HyperBall(G, (int) (Math.log(k_value)/Math.log(2)), 0);
    hyperball.init();

    for (int i = 1; i <= maxIter; i++) {
      hyperball.iterate();
      int modified = hyperball.modified();
      System.out.println("iteraton " + i + ": number of counters modified " + modified);
      if (modified == 0) {
        break;
      }
    }

    double effectiveD = NeighbourhoodFunction.effectiveDiameter(.9, hyperball.neighbourhoodFunction.toDoubleArray());
    hyperball.close();

    long t2 = System.currentTimeMillis();

    System.out.println("The effective diameter = " + effectiveD);
    System.out.println("Time elapsed (sec) for diameter computation = " + (t2 - t1)/1000.0);
  }
}
