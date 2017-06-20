import java.io.*;
import java.io.File;
import java.util.*;
import java.util.stream.*;

class Tuple {
  public double prob;
  public boolean label;

  public Tuple(double prob, boolean label) {
    this.prob = prob;
    this.label = label;
  }
}

public class Main {
  public static List<int[]> edgeSources;
  public static List<int[]> edgeTargets;
  public static List< List<Integer> > edgeTrainTest;  // [0,tp*2E): training; [tp*2E,2E): test
  public static List< Set<Integer> > nodes;
  public static List<int[]> freqs;		      // degree of each node, for average direction calculation

  public static Map<String, Integer> map = new HashMap<String, Integer>();
  public static Map<Integer, String> invMap = new HashMap<Integer, String>();

  /* system settings */
  public static int N;
  public static final int K = 4;
  public static final int NUM_OUTER_ITER = 2;	// number of outer iterations
  public static final int windowSize = 12;	// every time train on X consecutive timestamps (months)
  public static double scale = 2.0;		// range of initializations
  public static final double v = 0.03;		// speed
  public static double tp = 0.9;		// training portion (<1)
  public static double lr = 0.001;		// learing rate
  public static Random rand = new Random(0);

  /* default settings */
  public static double sigma = 1;
  public static double eps = 0.8;
  public static int MAX_ITER = 5000000;
  public static boolean VERBOSE = false;
  public static String prefix = "../../data/cosponsor/month_based/res/";
  public static int tStart = 120;
  public static int T = 100;			// time span: [tStart, tStart+T)

  /* model parameters */
  public static List<double[][]> xs = new ArrayList<double[][]>(T);	  // x: latent attribute (N*3: 2D position + 1D bias + 1D direction)


  /* pick a random number (id) from multiple bins */
  public static void linkSampler(int[] numEdges, int id, int[] res) {
    int bin = 0, edge = 0;
    try {
      while (edge + (int)(numEdges[bin] * tp) < id) {
	edge += (int)(numEdges[bin] * tp);
	bin++;
      }
      while (numEdges[bin] == 0) bin++;
    } catch (java.lang.ArrayIndexOutOfBoundsException e) {
      int sum = 0;
      for (int i: numEdges) {
	System.out.printf("%d\t", (int)(i*tp));
	sum += (int)(i*tp);
      }
      System.out.printf("\nsum = %d\t, id = %d\t%d\t%d\n", sum, id, bin, edge);
    }
    res[0] = bin; res[1] = id - edge;
  }

  /* parameter update: bias b */
  public static void train_b(int option) {
    List<int[]> sources = new ArrayList<int[]>(), targets = new ArrayList<int[]>();
    List< List<Integer> > edges = new ArrayList< List<Integer> >();
    int[] numEdges = new int[T];
    int sumEdges = 0; 
    for (int _t = 0; _t < T; _t++) {
      sources.add(edgeSources.get(_t));
      targets.add(edgeTargets.get(_t));
      edges.add(edgeTrainTest.get(_t));
      numEdges[_t] = edgeSources.get(_t).length;
      sumEdges += (int)(numEdges[_t] * tp);
    }

    for (int counter = 0; counter < MAX_ITER*10; counter++) {
      int[] res = new int[2];
      linkSampler(numEdges, rand.nextInt(sumEdges), res);
      int t = res[0];
      double[][] x0 = xs.get(0);
      double[][] x = xs.get(t);
      if (option == 0) x = x0;

      /* everytime sample a random link from G^{0:T} (training) */
      int e = edges.get(res[0]).get(res[1]);
      int E = numEdges[res[0]] / 2;
      int i = sources.get(res[0])[e], j = targets.get(res[0])[e];

      double[] grad = new double[6];

      double ins = -( (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) + (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) ) 
	/ (2.0 * eps * eps * x0[i][2] * x0[j][2]);
      if (e < E) {
	grad[0] += -( x[i][0]-x[j][0] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[1] += -( x[i][1]-x[j][1] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[2] += -( x[j][0]-x[i][0] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[3] += -( x[j][1]-x[i][1] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[4] += -ins / x0[i][2];
	grad[5] += -ins / x0[j][2];
      } else {
	if (ins > -1e-9) continue;
	double pij = Math.exp(ins);
	grad[0] += pij / (1-pij) * ( x[i][0]-x[j][0] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[1] += pij / (1-pij) * ( x[i][1]-x[j][1] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[2] += pij / (1-pij) * ( x[j][0]-x[i][0] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[3] += pij / (1-pij) * ( x[j][1]-x[i][1] ) / (eps * eps * x0[i][2] * x0[j][2]);
	grad[4] += pij / (1-pij) * ins / x0[i][2];
	grad[5] += pij / (1-pij) * ins / x0[j][2];
      }

      double thislr = 1e-4 * (1.0 - 1.0 * counter / (MAX_ITER * 10));
      if (t != 0 && option != 0) {
	x[i][0] += thislr * grad[0]; x[i][1] += thislr * grad[1];
	x[j][0] += thislr * grad[2]; x[j][1] += thislr * grad[3];
      } else {
	x0[i][0] += thislr * grad[0]; x0[i][1] += thislr * grad[1];
	x0[j][0] += thislr * grad[2]; x0[j][1] += thislr * grad[3];
      }
      /* update parameter b^{t} with projected gd */
      x0[i][2] += thislr * grad[4]; x0[j][2] += thislr * grad[5];
      if (x0[i][2] < 1e-9) x0[i][2] = 1e-9;
      if (x0[j][2] < 1e-9) x0[j][2] = 1e-9;

      if (t != 0 && option != 0) xs.set(t, x);
      xs.set(0, x0);

      /*
      if ((counter+1) % MAX_ITER == 0) {
	if (t != 0 && option != 0) {
	  double new_obj = calcObj(0, T, x);
	  System.out.printf("[Training] counter = %d\tobj = %f\n", counter+1, new_obj);
	} else {
	  double new_obj = calcObj(0, T, x0);
	  System.out.printf("[Training] counter = %d\tobj = %f\n", counter+1, new_obj);
	}
      }
      */
    }

    if (option == 0) {
      double[][] x0 = xs.get(0);
      double[][] new_x = new double[N][K];
      for (int n = 0; n < N; n++) {
	/* update x_{i/j}^{t} according to v^{t} */
	new_x[n][0] = x0[n][0] + v * CosSineTable.getCos(x0[n][3]);
	new_x[n][1] = x0[n][1] + v * CosSineTable.getSine(x0[n][3]);
	new_x[n][2] = x0[n][2];
	new_x[n][3] = x0[n][3];
      }
      xs.set(0, x0);
      xs.set(1, new_x);
      System.out.print("\n\n");
      FileParser.output(x0, invMap, "./res/x_0.txt", new HashSet<Integer>());
    }

    return;
  }

 
  public static void train() {
    /* read, init data & parameters */
    N = FileParser.readCSVDict(prefix, tStart, T, map, invMap);
    for (int t = 0; t < T; t++) {
      String fileName = prefix + "/" + Integer.toString(tStart+t) + ".csv";

      int E = FileParser.countLines(fileName);
      int[] edgeSource = new int[2*E], edgeTarget = new int[2*E];
      Set<Integer> curNodes = new HashSet<Integer>();
      int[] freq = new int[N];
      double nsw = FileParser.readCSVGraph(fileName, N, E, map, edgeSource, edgeTarget, curNodes, freq);
      edgeSources.add(edgeSource); edgeTargets.add(edgeTarget);

      int curN = curNodes.size();
      nodes.add(curNodes);
      freqs.add(freq);

      List<Integer> allEdges = new ArrayList<Integer>(2*E);
      for (int e = 0; e < 2*E; e++) allEdges.add(e);
      Collections.shuffle(allEdges, new Random(t));
      edgeTrainTest.add(allEdges);

      double[][] x = new double[N][K];
      for (int n = 0; n < N; n++) {
	x[n][0] = scale * (rand.nextDouble() - 0.5);	    // coord 1
	x[n][1] = scale * (rand.nextDouble() - 0.5);	    // coord 2
	x[n][2] = 1 + 0.1 * rand.nextDouble();		    // bias b
	x[n][3] = 2 * Math.PI * (rand.nextDouble() - 0.5);  // direction
      }
      xs.add(x);
    }

    /* training: 
     * update 
     *	  theta_{t-1} (and x_{t} = x_{t-1} + v theta_{t-1}) 
     * given 
     *	  G_{t}, theta_{t-2}, x_{t-1} (index starts from 0)
     *
     * (objective evaluated on G_{t})
     */

    for (int outer_iter = 0; outer_iter < NUM_OUTER_ITER; outer_iter++) {
      if (VERBOSE) System.out.printf("[Training] Outer iteration = %d\n", outer_iter);
      train_b(outer_iter);

      for (int t = 1; t < T-windowSize; t++) {
	Set<Integer> curNodes = getUnionNodes(t-1, t+windowSize-1);
	double[][] x = xs.get(t-1);	// x^{t-1}
	double[][] new_x = xs.get(t);	// x^{t}

	if (t > 1) {
	  double[][] pre_x = xs.get(t-2);	// x^{t-2}
	  double[] averDir = averDirection(t-2, t-2+windowSize);
	  for (int n = 0; n < N; n++) {
	    for (int k = 0; k < 4; k++) {
	      new_x[n][k] = x[n][k];
	    }
	  }
	  for (int n: curNodes) {
	    x[n][3] = averDir[n];
	    new_x[n][0] = x[n][0] + v * CosSineTable.getCos(x[n][3]);
	    new_x[n][1] = x[n][1] + v * CosSineTable.getSine(x[n][3]);
	  }
	}

	double old_obj = calcObj(t, t+windowSize, new_x);
	if (VERBOSE) System.out.printf("t from %d to %d:\tobj (init) = %f\n", t, t+windowSize, old_obj);

	List<int[]> sources = new ArrayList<int[]>(), targets = new ArrayList<int[]>();
	List< List<Integer> > edges = new ArrayList< List<Integer> >();
	int[] numEdges = new int[windowSize];
	int sumEdges = 0; 
	for (int _t = 0; _t < windowSize; _t++) {
	  sources.add(edgeSources.get(t+_t));
	  targets.add(edgeTargets.get(t+_t));
	  edges.add(edgeTrainTest.get(t+_t));
	  numEdges[_t] = edgeSources.get(t+_t).length;
	  sumEdges += (int)(numEdges[_t] * tp);
	}

	for (int counter = 0; counter < MAX_ITER; counter++) {
	  int[] res = new int[2];
	  linkSampler(numEdges, rand.nextInt(sumEdges), res);
	  /* everytime sample a random link (training) */
	  int e = edges.get(res[0]).get(res[1]);
	  int E = numEdges[res[0]] / 2;
	  int i = sources.get(res[0])[e], j = targets.get(res[0])[e];
	  double grad_theta_i = 0, grad_theta_j = 0;	  // gradient for theta (3)

	  /* first term: p(G^{t} | x^{t}) = p(G^{t} | x^{t-1}, theta^{t-1}) */
	  double ins = -( (new_x[i][0]-new_x[j][0]) * (new_x[i][0]-new_x[j][0]) 
	      + (new_x[i][1]-new_x[j][1]) * (new_x[i][1]-new_x[j][1]) ) 
	    / (2.0 * eps * eps * new_x[i][2] * new_x[j][2]);
	  if (e < E) {
	    grad_theta_i += ( (new_x[i][0]-new_x[j][0]) * CosSineTable.getSine(x[i][3]) 
		- (new_x[i][1]-new_x[j][1]) * CosSineTable.getCos(x[i][3]) )
	      * v / (eps * eps * new_x[i][2] * new_x[j][2]);
	    grad_theta_j += ( (new_x[j][0]-new_x[i][0]) * CosSineTable.getSine(x[j][3]) 
		- (new_x[j][1]-new_x[i][1]) * CosSineTable.getCos(x[i][3]) )
	      * v / (eps * eps * new_x[i][2] * new_x[j][2]);
	  } else {
	    if (ins == 0) continue;
	    double pij = Math.exp(ins);
	    grad_theta_i -= ( (new_x[i][0]-new_x[j][0])* CosSineTable.getSine(x[i][3]) 
		- (new_x[i][1]-new_x[j][1]) * CosSineTable.getCos(x[i][3]) ) 
	      * v / (eps * eps * new_x[i][2] * new_x[j][2]) * pij / (1-pij);
	    grad_theta_j -= ( (new_x[j][0]-new_x[i][0]) * CosSineTable.getSine(x[j][3]) 
		- (new_x[j][1]-new_x[i][1]) * CosSineTable.getCos(x[i][3]) )
	      * v / (eps * eps * new_x[i][2] * new_x[j][2]) * pij / (1-pij);
	  }

	  double thislr = lr * (1.0 - 1.0 * counter / MAX_ITER);
	  /* update parameter v^{t-1} */
	  x[i][3] += thislr * grad_theta_i; 
	  x[j][3] += thislr * grad_theta_j;

	  /* second term: p(theta^{t-1} | theta^{t-2}) - only valid when t > 1 */
	  if (counter%(MAX_ITER/10) == 0 && t > 1) {
	    double[] averDir = averDirection(t-2, t-2+windowSize);
	    for (int n: curNodes) {	// only those n which appear in both G^{t-2} and G^{t-1} are considered
	      Set<Integer> preNodes = getUnionNodes(t-2, t+windowSize-2);
	      if (preNodes.contains(n)) {
		double ad = Util.angleDistance(x[n][3], averDir[n]);
		/* update parameter v^{t-1}, and x_{i/j}^{t} according to v^{t} */
		x[n][3] -= thislr * ad / (sigma*sigma);
		new_x[n][0] = x[n][0] + v * CosSineTable.getCos(x[n][3]);
		new_x[n][1] = x[n][1] + v * CosSineTable.getSine(x[n][3]);
	      }
	    }
	  }

	  /* update x_{i/j}^{t} according to v^{t} */
	  new_x[i][0] = x[i][0] + v * CosSineTable.getCos(x[i][3]);	  
	  new_x[i][1] = x[i][1] + v * CosSineTable.getSine(x[i][3]);
	  new_x[j][0] = x[j][0] + v * CosSineTable.getCos(x[j][3]);
	  new_x[j][1] = x[j][1] + v * CosSineTable.getSine(x[j][3]);

	  /*
	  if ((counter+1) % (MAX_ITER/10) == 0) {
	    double new_obj = calcObj(t, t+windowSize, new_x);
	    System.out.printf("[Training] counter = %d\tobj = %f\n", counter+1, new_obj);
	    double rate = (new_obj-old_obj) / Math.abs(old_obj);
	    //if (rate < 1e-6) break;
	    old_obj = new_obj;
	  }
	  */
	}

	/* update v^{t-1} (x[3] has changed) and x^{t} (new_x[0:2] has changed) */
	xs.set(t-1, x);
	xs.set(t, new_x);
	File _f = new File("./res_" + outer_iter);
	_f.mkdir();
	FileParser.output(x, invMap, "./res_" + outer_iter + "/x_" + (t-1) + ".txt", curNodes);
      }
    }	// end outer_iter
  }

  /* calcObj:
   *	graph: [t1, t2)
   *	objective = log p( G^{ t1 : (t2-1) } | x^{t1} )	      + log p(theta^{t1} | theta^{ t1-1 }) 
   */
  public static double calcObj(int t1, int t2, double[][] x) {
    double res = 0;
    int numEdges = 0;
    for (int t = t1; t < t2; t++) {
      int[] edgeSource = edgeSources.get(t);
      int[] edgeTarget = edgeTargets.get(t);
      List<Integer> allEdges = edgeTrainTest.get(t);
      int E = edgeSource.length / 2;
      numEdges += (int)(tp*2*E);

      for (int _e = 0; _e < tp*2*E; _e++) {
	int e = allEdges.get(_e);
	int i = edgeSource[e], j = edgeTarget[e];
	double ins = -( (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) + (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) )
	  / (2.0 * eps * eps * x[i][2] * x[j][2]);
	if (e < E) {
	  res += ins;
	} else if (ins < -1e-9) {
	  res += Math.log(1-Math.exp(ins));
	}
      }
    }

    if (t1 > 0) {
      Set<Integer> curNodes = getUnionNodes(t1, t2);
      double[] averDir = averDirection(t1-1, t2-1);
      Set<Integer> preNodes = getUnionNodes(t1-1, t2-2);
      for (int n: curNodes) {	// only those n which appear in both G^{t-2} and G^{t-1} are considered
	if (preNodes.contains(n)) {
	  double ad = Util.angleDistance(x[n][3], averDir[n]);
	  res -= 0.5 * (ad/sigma) * (ad/sigma);
	}
      }
    }

    return res/numEdges;
  }


  /* average direction of neighbors at time t */
  public static double[] averDirection(int t1, int t2) {
    double[][] pos = new double[N][2];
    for (int t = t1; t < t2; t++) {
      int[] edgeS = edgeSources.get(t), edgeT = edgeTargets.get(t);
      int E = edgeS.length / 2;
      double[][] x = xs.get(t);
      List<Integer> allEdges = edgeTrainTest.get(t);
      for (int _e = 0; _e < tp*2*E; _e++) {
	int e = allEdges.get(_e);
	if (e < E) {
	  int i = edgeS[e], j = edgeT[e];
	  pos[i][0] += CosSineTable.getCos(x[j][3]);
	  pos[i][1] += CosSineTable.getSine(x[j][3]);
	  pos[j][0] += CosSineTable.getCos(x[i][3]);
	  pos[j][1] += CosSineTable.getSine(x[i][3]);
	}
      }
    }

    double[] res = new double[N];
    for (int n = 0; n < N; n++) {
      res[n] = Math.atan2(pos[n][1], pos[n][0]);    // range of Math.atan2 is [-Pi, Pi]
    }
    return res;
  }


  /* get active nodes from [t1, t2) */
  public static Set<Integer> getUnionNodes(int t1, int t2) {
    Set<Integer> res = new HashSet<Integer>();
    for (int t = t1; t < t2; t++) {
      Set<Integer> node_t = nodes.get(t);
      for (int s: node_t) res.add(s);
    }
    return res;
  }


  public static void main(String[] args) {
    try {
      ArgumentParser.parse(args);
    } catch (NumberFormatException e) {
      //e.printStackTrace();
      System.out.println("Illegal arguments.");
      ArgumentParser.help();
      System.exit(0);
    }

    edgeSources = new ArrayList<int[]>(T);
    edgeTargets = new ArrayList<int[]>(T);
    edgeTrainTest = new ArrayList< List<Integer> >(T);
    nodes = new ArrayList< Set<Integer> >(T);
    freqs = new ArrayList<int[]>(T);

    train();
  }

}


