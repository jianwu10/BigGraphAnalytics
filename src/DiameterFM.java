/*
 * Diameter and effective diameter estimation using the Flajolet-Martin algorithm.
 * Reference: https://epubs.siam.org/doi/pdf/10.1137/1.9781611972801.48
 */

import java.util.ArrayList;
import java.util.Arrays;
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;

public class DiameterFM {
  ImmutableGraph G;
  int n; // number of nodes
  int K; // number of bitstrings: i.e. 32
  int MAX_ITER; // max number of iterations: ~ 256
  int diameter;
  double effectiveD;
  int[] bitstrings_curr; // current bitstring array for the nodes
  int[] bitstrings_prev; // bitstring array in the previous iteration

  public DiameterFM(String basename, int K, int MAX_ITER) throws Exception {
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

  public void computeDiamter()  {
    int hMax = MAX_ITER;
    AtomicInteger changed = new AtomicInteger(0);
    ArrayList<Double> N = new ArrayList<>();
    // h = 0, N[0] = n
    N.add((double) n);
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
      // update bitstrings_curr using bit OR operation
      Arrays.stream(nodeIndex)
            .parallel()
            .boxed()
            .mapToInt(u -> {
              ImmutableGraph H = G.copy();
              int[] neighbors = H.successorArray(u);
              for (int i = 0; i < K; i++) {
            	  for (int neighbor : neighbors) {
            	    bitstrings_curr[u*K+i] = bitstrings_curr[u*K+i] | bitstrings_prev[neighbor*K+i];
            	  }
            	  if (bitstrings_curr[u*K+i] != bitstrings_prev[u*K+i]) {
                  changed.getAndIncrement();
            	  }
            	}
            	return 0;
            }).sum();
      // append the sum of the neighborhood function output to N
      double sum = Arrays.stream(nodeIndex)
                        .parallel()
                        .boxed()
                        .mapToDouble(u -> {
                          double average = 0;
                          for (int i = 0; i < K; i++) {
                            average += firstZeroPos(bitstrings_curr[u*K+i]);
                          }
                          average /= K;
                          return Math.pow(2, average) / 0.77351;
                        }).sum();
      N.add(sum);
      System.out.println("number of changed bitstrings = " + changed);
      if (changed.get() == 0) {
        hMax = h;
        break;
      }
    }
    // Now we have hmax which is the diameter.
    diameter = hMax;
    // Compute the effective diameter.
    // The effective diameter is the smallest h where N(h) = 0.9 N(hMax).
    int index_N;
    for (index_N = 1; index_N < N.size(); index_N++) {
      double ninetyPercent = N.get(N.size()-1) * 0.9;
      if (N.get(index_N) >= ninetyPercent) {
        double fraction = (ninetyPercent - N.get(index_N-1))/(N.get(index_N) - N.get(index_N-1));
        effectiveD = index_N  - 1 + fraction;
        break;
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
    DiameterFM DFM = new DiameterFM("data/"+basename+"/"+basename, k_value, maxIter);
    DFM.computeDiamter();
    long t2 = System.currentTimeMillis();
    System.out.println("The diameter = " + DFM.diameter);
    System.out.println("The effective diameter = " + DFM.effectiveD);
    System.out.println("Time elapsed (sec) for diameter computation = " + (t2 - t1)/1000.0);
  }
}
