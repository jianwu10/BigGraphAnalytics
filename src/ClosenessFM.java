/*
 * Closeness centrality estimation using the Flajolet-Martin algorithm.
 * Reference: http://igraph.org/r/doc/closeness.html
 */

import java.util.ArrayList;
import java.util.Arrays;
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.File;
import java.io.PrintStream;

public class ClosenessFM {
  ImmutableGraph G;
  int n; // number of nodes
  int K; // number of bitstrings: i.e. 32
  int MAX_ITER; // max number of iterations: ~ 256
  int[] bitstrings_curr; // current bitstring array for the nodes
  int[] bitstrings_prev; // bitstring array in the previous iteration
  double[] N_curr; // neighborhodd function value in current hop for each node
  double[] sumOfDist; // sum of distances for each node
  double[] closeness;

  public ClosenessFM(String basename, int K, int MAX_ITER) throws Exception {
    System.out.println("Starting initialization.");
    G = ImmutableGraph.loadMapped(basename);
    n = G.numNodes();
    this.K = K;
    this.MAX_ITER = MAX_ITER;
    bitstrings_curr = new int[n*this.K];
    bitstrings_prev = new int[n*this.K];
    int[] nodeIndex = new int[n];
    Arrays.parallelSetAll(nodeIndex, i -> i);
    Arrays.stream(nodeIndex)
          .parallel()
          .boxed()
          .mapToInt(u -> {
            for (int i = 0; i < this.K; i++) {
              bitstrings_curr[u*this.K+i] = 1 << bitIndex();
            }
            return 0;
          }).sum();
    System.out.println("Finished.");
    N_curr = new double[n];
    closeness = new double[n];
    Arrays.parallelSetAll(N_curr, i -> 1.0);
    sumOfDist = new double[n];
  }

  public int bitIndex() {
    int i = 0;
    double randomNum = Math.random();
    double threshold = 0;
    for (i = 0; i < 32; i++) {
      threshold += Math.pow(2, -1*i-1);
      if (randomNum < threshold) {
        break;
      }
    }
    return i;
  }

  // find the rightmost "0" position
  public int firstZeroPos(int bs) {
    int i = 0;
    for (i = 0; i < 32; i++) {
      int flag = bs & (1 << i);
      if (flag == 0) {
        return i;
      }
    }
    return i;
  }

  public void computeCloseness()  {
    AtomicInteger changed = new AtomicInteger(0);
    int[] nodeIndex = new int[n];
    Arrays.parallelSetAll(nodeIndex, i -> i);
    for (int h = 1; h <= MAX_ITER; h++) {
      System.out.println("Starting iteration " + h);
      // reset changed for this iteration
      changed.set(0);
      // update bitstrings_prev
      for (int i = 0; i < bitstrings_prev.length; i++) {
        bitstrings_prev[i] = bitstrings_curr[i];
      }
      final int h_curr = h;
      Arrays.stream(nodeIndex)
            .parallel()
            .boxed()
            .mapToInt(u -> {
              ImmutableGraph H = G.copy();
              int[] neighbors = H.successorArray(u);
              // update N_prev
              double N_prev = N_curr[u];
              double average = 0;
              for (int i = 0; i < K; i++) {
            	  for (int neighbor : neighbors) {
                  // update bitstrings_curr using bit OR operation
            	    bitstrings_curr[u*K+i] = bitstrings_curr[u*K+i] | bitstrings_prev[neighbor*K+i];
            	  }
                average += firstZeroPos(bitstrings_curr[u*K+i]);
            	  if (bitstrings_curr[u*K+i] != bitstrings_prev[u*K+i]) {
                  changed.getAndIncrement();
            	  }
            	}
              average /= K;
              // update N_curr
              N_curr[u] = Math.pow(2, average) / 0.77351;
              // update sumOfDist
              if (N_curr[u] > N_prev) {
                sumOfDist[u] += h_curr * (N_curr[u] - N_prev);
              }
            	return 0;
            }).sum();
      System.out.println("number of changed bitstrings = " + changed);
      if (changed.get() == 0) {
        break;
      }
    }
    // Calcuate closeness based on sumOfDist
    // Normalize closeness by n-1.
    // For those unreachable nodes, use n as the distance.
    for (int i = 0; i < n; i++) {
      double sum = sumOfDist[i];
      if (N_curr[i] < n) {
        // There are unreachable nodes for i
        sum += n * (n - N_curr[i]);
        closeness[i] = (n - 1) / sum;
      } else {
        closeness[i] = (n - 1) / sum;
      }
    }
  }

  public static void main(String[] args) throws Exception {
    System.out.println("Please enter graph's base name");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();

    System.out.println("Please enter the K value");
    int k_value = sc.nextInt();

    System.out.println("Please enter the max iteration number");
    int maxIter = sc.nextInt();

    long t1 = System.currentTimeMillis();
    ClosenessFM CFM = new ClosenessFM("data/"+basename+"/"+basename, k_value, maxIter);
    CFM.computeCloseness();
    long t2 = System.currentTimeMillis();
    System.out.println("Time elapsed (sec) for closeness computation = " + (t2 - t1)/1000.0);

    // define path strings
    String inputPath = "data/";
    String outputPath = "/output/closenessCentrality/";
    PrintStream ps = new PrintStream(new File(inputPath+basename+"/"+outputPath+basename+"_closenessCentrality_FM_.txt"));
    System.out.println("writing results to file...");
    for (int i = 0; i < CFM.closeness.length; i++) {
      ps.println(i + "\t" + CFM.closeness[i]);
    }

    System.out.println("Sorting... be patient");
    Runtime.getRuntime().exec("sort -g -r -k2,2 -k1,1n "+inputPath+basename+"/"+outputPath+basename+"_closenessCentrality_FM_.txt -o "+inputPath+basename+"/"+outputPath+basename+"_closenessCentrality_FM_sorted.txt").waitFor();
  }
}
