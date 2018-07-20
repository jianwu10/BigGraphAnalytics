/*
*  H-index based k-truss decomposition https://arxiv.org/pdf/1704.00386.pdf
*  Synchronous
*/
import java.util.List;
import java.io.File;
import java.io.PrintStream;
import java.io.FileOutputStream;
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicIntegerArray;

public class KtrussSync {
  final ImmutableGraph G; // WebGraph object
  final int n; // number of nodes
  final int m; // number of edges
  final int[] edgeTail;
  final int[] offset;  // edge offset for edgeHead array
  final int[] degree;
  AtomicIntegerArray supportCurr; // support for each edge for current iteration
  int[] supportPrev; // support for each edge for previous iteration
  BitSet updated; // to record is an edge is updated or not

  public KtrussSync(String basename) throws Exception {
    G = ImmutableGraph.loadMapped(basename);
    n = G.numNodes();
    degree = new int[n];
    offset = new int[n+1];
    // Since we have preprocessed the graph, there are no self-loops anymore.
    // m is the number of edges in the undirected graph.
    m = (int) (G.numArcs()/2);
    edgeTail = new int[m];
    // Populate edgeHead and edgeTail such that u < v.
    System.out.println("Start populating edgeHead and edgeTail.");
    int startPos = 0;
    for (int u = 0; u < n; u++) {
      offset[u] = startPos;
      // Get the neighbors of u
      degree[u] = G.outdegree(u);
      int[] u_successor = G.successorArray(u);
      int count = 0;
      for (int j = 0; j < u_successor.length; j++) {
        int v = u_successor[j];
        if (u <  v) {
          edgeTail[startPos+count] = v;
          count++;
        }
      }
      startPos += count;
    }
    offset[n] = offset[n-1];
    System.out.println("Finished.");
    System.out.println("Total number of edges = " + m);
    supportCurr = new AtomicIntegerArray(m);
    supportPrev = new int[m];
    updated = new BitSet(m);
  }

  public void supportCompute() {
    // Step 1: compute support for each edge.
    System.out.println();
    System.out.println("Starting step 1: support computing for each edge.");
    long startTime1 = System.currentTimeMillis();
    // Since m is huge, we divide the edge list into segments.
    // We parallelize each segment at a time.
    int numOfSeg = 10000;
    int segLength = m/numOfSeg;
    int residual = m%numOfSeg;
    int[] index = new int[segLength];
    for (int segCount = 0; segCount < numOfSeg; segCount++) {
      int start = segCount * segLength;
      // Each index is the index in edgeHead (or edgeTail) array.
      Arrays.parallelSetAll(index, (i) -> (start+i));
      Arrays.stream(index)
            .parallel()
            .boxed()
            .mapToInt(e -> {
              int u = findHead(e);
              int v = edgeTail[e];
              // We need to identify triangles containing edge uv.
              // We define a triangle to be uvw such that u < v < w.
              // We just need to find all w.
              // Do the intersection of w1 and w2
              int w1 = e + 1;
              int w2 = offset[v];
              while (w1 < offset[u+1] && w2 < offset[v+1]) {
                if (edgeTail[w1] < edgeTail[w2]) {
                  w1++;
                } else if (edgeTail[w1] > edgeTail[w2]) {
                  w2++;
                } else {
                  // edgeTail[w1] = edgeTail[w2]
                  // Now we find a triangle.
                  // We add 1 to edge uv, vw, and uw's support
                  supportCurr.getAndIncrement(e);
                  supportCurr.getAndIncrement(w1);
                  supportCurr.getAndIncrement(w2);
                  w1++;
                  w2++;
                }
              }
              return 0;
      }).sum();
      // System.out.println("Calculated support for " + (segCount+1)*100.0/numOfSeg + "% edges.");
    }
    // Calculate support for the residual edges
    for (int k = 0; k < residual; k++) {
      int start = m - residual;
      int e = start + k;
      int u = findHead(e);
      int v = edgeTail[e];
      // We need to identify triangles containing edge uv.
      // We define a triangle to be uvw such that u < v < w.
      // We just need to find all w.
      // Do the intersection of w1 and w2
      int w1 = e + 1;
      int w2 = offset[v];
      while (w1 < offset[u+1] && w2 < offset[v+1]) {
        if (edgeTail[w1] < edgeTail[w2]) {
          w1++;
        } else if (edgeTail[w1] > edgeTail[w2]) {
          w2++;
        } else {
          // edgeTail[w1] = edgeTail[w2]
          // Now we find a triangle.
          // We add 1 to edge uv, vw, and uw's support
          supportCurr.getAndIncrement(e);
          supportCurr.getAndIncrement(w1);
          supportCurr.getAndIncrement(w2);
          w1++;
          w2++;
        }
      }
    }
    System.out.println("maximum support: " + findMaxSupport(supportCurr));
    long endTime1 = System.currentTimeMillis();
    System.out.println("Time elapsed (sec) for support computing = " + (endTime1 - startTime1)/1000.0);
  }

  public void ktrussCompute() {
    // Step 2: k-truss computation (synchronous)
    System.out.println();
    System.out.println("Starting step 2: k-truss computation (synchronous).");
    long startTime2 = System.currentTimeMillis();
    int iter = 0;
    // Set one bit to start the loop
    updated.set(0);
    // While the updated is ture
    while (!updated.isEmpty()) {
      iter++;
      System.out.print("iteration " + iter + ": ");
      // Set updated to false
      updated.clear();
      // Pass supportCurr's values to supportPre
      for (int i = 0; i < m; i++) {
        supportPrev[i] = supportCurr.get(i);
      }
      // For each edge, we calculate the h-index in parallel.
      // We just parallelize each segment to avoid out-of-memory exception.
      // Calculate the h-index for each edge in each segment.
      int numOfSeg = 10000;
      int segLength = m/numOfSeg;
      int residual = m%numOfSeg;
      int[] index = new int[segLength];
      for (int segCount = 0; segCount < numOfSeg; segCount++) {
        int start = segCount * segLength;
        // Each index is the index in edgeHead (or edgeTail) array.
        Arrays.parallelSetAll(index, (i) -> (start+i));
        Arrays.stream(index)
              .boxed()
              .parallel()
              .mapToInt(e -> {
                int u = findHead(e);
                int v = edgeTail[e];
                // Make sure that degree(u) < degree(v)
                if (degree[u] > degree[v]) {
                  // Swap u and v
                  int temp = u;
                  u = v;
                  v = temp;
                }
                // Get the neighbors of u
                ImmutableGraph H = G.copy();
                int[] u_neighbor = H.successorArray(u);
                ArrayList<Integer> L = new ArrayList<>(degree[u]);
                // For each neighbor of u
                for (int w : u_neighbor) {
                  // Check if <v, w> is an edge or not
                  // If <v,w> is an edge, get its id.
                  int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u,w)], offset[Math.min(u,w)+1], Math.max(u,w));
                  int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v,w)], offset[Math.min(v,w)+1], Math.max(v,w));
                  if (id_vw > 0) {
                    int rho = Math.min(supportPrev[id_uw], supportPrev[id_vw]);
                    L.add(rho);
                  }
                }
                int hIndex = getHIndex(L);
                if (hIndex != supportPrev[e]) {
                  updated.set(e);
                }
                supportCurr.getAndSet(e, hIndex);
                return 0;
              }).sum();
      }
      // Calculate h-index for the residual edges
      for (int i = 0; i < residual; i++) {
        int start = m - residual;
        int e = start + i;
        int u = findHead(e);
        int v = edgeTail[e];
        // Make sure that degree(u) < degree(v)
        if (degree[u] > degree[v]) {
          // Swap u and v
          int temp = u;
          u = v;
          v = temp;
        }
        // Get the neighbors of u
        int[] u_neighbor = G.successorArray(u);
        ArrayList<Integer> L = new ArrayList<>(degree[u]);
        // For each neighbor of u
        for (int w : u_neighbor) {
          // Check if <v, w> is an edge or not
          // If <v,w> is an edge, get its id.
          int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u,w)], offset[Math.min(u,w)+1], Math.max(u,w));
          int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v,w)], offset[Math.min(v,w)+1], Math.max(v,w));
          if (id_vw > 0) {
            int rho = Math.min(supportPrev[id_uw], supportPrev[id_vw]);
            L.add(rho);
          }
        }
        int hIndex = getHIndex(L);
        if (hIndex != supportPrev[e]) {
          updated.set(e);
        }
        supportCurr.getAndSet(e, hIndex);
      }
      // Now we have finished one iteration to calculate h-index for each eage.
      System.out.println("Number of edges updated = " + updated.cardinality());
    }

    long endTime2 = System.currentTimeMillis();
    System.out.println("Time elapsed (sec) for k-truss computation = " + (endTime2 - startTime2)/1000.0);

    // Find the maximum truss
    // kmax is the maximum truss.
    // Or we can find the maximum support in the current support array.
    // Then max truss = max support + 2
    int maxTruss = findMaxSupport(supportCurr) + 2;
    System.out.println("Max truss = " + maxTruss);
  }

  // Method to calculate h-index for an array
  private int getHIndex(ArrayList<Integer> L) {
    // O(n) time. Using bucket.
    if (L == null || L.size() == 0) {
      return 0;
    }
    int n = L.size();
    int[] bucket = new int[n+1];
 		// Fill the bucket
 		for (int i = 0; i < n; i++) {
      int bucketSlot = L.get(i);
 			if (L.get(i) <= n) {
 			  bucket[bucketSlot]++;
 			} else {
 				bucket[n]++;
 			}
 		}
 		//Find the h-index
 		int sum = 0;
 		for (int h = n; h >= 0; h--) {
 			if (sum + bucket[h] >= h) {
 				return h;
 			}
 			sum += bucket[h];
 		}
    return 0;
  }

  private int findMaxSupport(AtomicIntegerArray support) {
    int max = 0;
    for (int i = 0; i < support.length(); i++) {
      int temp = support.get(i);
      if (temp > max) {
        max = temp;
      }
    }
    return max;
  }

  private int findHead(int target) {
    // find index such that offset[index] <= target
    // and target < offset[index+1]
    int start = 0;
    int end = offset.length - 1;
    while (start + 1 < end) {
      int mid = start + (end - start)/2;
      if (offset[mid] <= target) {
        start = mid;
      } else {
        end = mid;
      }
    }
    if (offset[end] <= target) {
      return end;
    }
    if (offset[start] <= target) {
      return start;
    }
    return 0;
  }

  public static void main(String[] args) throws Exception{
    System.out.println("Please enter graph's base name");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();

    // PrintStream out = new PrintStream(new FileOutputStream("data/"+basename+"/"+"output/ktruss/"+basename+"_ktruss_synchronous_console.txt"));
    // System.setOut(out);

    KtrussSync kt = new KtrussSync("data/"+basename+"/"+basename+"-noself-sym");
    kt.supportCompute();
    kt.ktrussCompute();

    // System.out.println("Writing results to file. Be patient...");
    // PrintStream ps = new PrintStream(new File("data/"+basename+"/"+"output/ktruss/"+basename+"_ktruss_synchronous.txt"));
    // for (int i = 0; i < kt.m; i++) {
    //   ps.println(kt.supportCurr.get(i)+2);
    // }
  }
}
