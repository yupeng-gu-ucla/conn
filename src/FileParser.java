import java.io.*;
import java.util.*;
import java.lang.*;

public class FileParser {
  public static final String delim = ",";
  public static final int NEG_TABLE_SIZE = 1000000;

  public static int countLines(String filename) {
    int count = 0; 
    boolean empty = true; 
    try (InputStream is = new BufferedInputStream(new FileInputStream(filename))) {
      byte[] c = new byte[1024]; 
      int readChars = 0; 
      while ((readChars = is.read(c)) != -1) { 
	empty = false; 
	for (int i = 0; i < readChars; ++i) { 
	  if (c[i] == '\n') {
	    ++count;
	  }
	}
      }
    } catch (IOException e) {}
    return (count == 0 && !empty) ? 1 : count;
  }

  /**
   * readCSVDict: 
   *	build a map (original ID: new ID in memory) and an inverted map
   */
  public static int readCSVDict(String prefix, int tStart, int T, Map<String, Integer> map, Map<Integer, String> invMap) {
    int newID = 0;
    for (int t = 0; t < T; t++) {
      String filename = prefix + "/" + Integer.toString(tStart+t) + ".csv";

      try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
	String currentLine;
	int lineCount = 0;
	while ((currentLine = br.readLine()) != null) {
	  //if (lineCount % 1000000 == 0) System.out.println(lineCount);
	  lineCount += 1;

	  String[] tokens = currentLine.split(delim);
	  String x = tokens[0];
	  String y = tokens[1];
	  if (!map.containsKey(x)) {
	    map.put(x, newID); newID += 1;
	  }
	  if (!map.containsKey(y)) {
	    map.put(y, newID); newID += 1;
	  }
	}
      } catch (IOException _e) {
	//System.out.printf("[Info] File %s does not exist (readCSVDict).\n", filename);
      } 

      // set inverse map
      for (Map.Entry<String, Integer> e: map.entrySet()) {
	invMap.put(e.getValue(), e.getKey());
      }
    }
    //System.out.printf("[Info] ID count = %d, map size = %d\n", newID, map.size());

    try (Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("ID_map.txt"), "utf-8"))) {
      for (Map.Entry<String, Integer> e: map.entrySet()) {
	writer.write(String.format("%s\t%d\n", e.getKey(), e.getValue()));
      }
    } catch (IOException e) {
      System.out.printf("[Info] Mapping output error\n");
    }

    return map.size();
  }

  /** readCSVGraph:
   *	read network in .csv format
   *	A: inverse of out-degree
   *	label: whether the file contains existing links or non-existing links
   *  @return:
   *	negative sample weight
   */
  public static double 
  readCSVGraph(String fileDir, int N, int E, Map<String, Integer> map, 
      int[] edgeSource, int[] edgeTarget, 
      Set<Integer> curNodes, int[] freq
  ) {
    Set<Long> pairs = new HashSet<Long>();
    try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
      String currentLine;
      int e = 0;
      while ((currentLine = br.readLine()) != null) {
	String[] tokens = currentLine.split(delim);
	int x = map.get(tokens[0]);
	int y = map.get(tokens[1]);
	curNodes.add(x); curNodes.add(y);
	edgeSource[e] = x; edgeTarget[e] = y;
	freq[x]++; 
	freq[y]++;	// undirected
	e++;

	long pairKey = x * N + y;
	pairs.add(pairKey);
      }

      // negative samples
      int[] negTable = new int[NEG_TABLE_SIZE];

      int localN = curNodes.size();
      double sum = 0, por = 0, cur_sum = 0;
      int[] nodes = new int[localN];
      int i = 0;
      for (int n: curNodes) {
	nodes[i] = n;
	i++;
      }
      for (i = 0; i < localN; i++) {
	int n = nodes[i];
	sum += Math.pow(freq[n], 0.75);
      }
      i = 0;
      for (int k = 0; k < NEG_TABLE_SIZE; k++) {
	if (1.0 * (k+1) / NEG_TABLE_SIZE > por) {
	  int n = nodes[i];
	  cur_sum += Math.pow(freq[n], 0.75);
	  por = cur_sum / sum;
	  i++;
	}
	negTable[k] = nodes[i-1];
      }

      Random rand = new Random(0);
      for (e = 0; e < E; e++) {
	int y = edgeTarget[e];
	//int x = edgeSource[e];
	while (true) {
	  int x = nodes[rand.nextInt(localN)];

	  //int y = nodes[rand.nextInt(localN)];
	  ////////int y = negTable[rand.nextInt(NEG_TABLE_SIZE)];
	  //if (!pairs.contains(x * N + y)) {
	  if (true) {
	    edgeSource[E+e] = x; edgeTarget[E+e] = y;
	    break;
	  }
	}
      }
    } catch (IOException _e) {
      //System.out.printf("[Info] File %s does not exist (readCSVgraph).\n", fileDir);
      //_e.printStackTrace();
    }

    if (pairs.size() == 0) return 0;
    else return 1.0*N*(N-1)/pairs.size()-1;
  }

  public static void
  output(double[][] arr, Map<Integer, String> inv_map, String fileDir, Set<Integer> curNodes) {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
      for (int i = 0; i < arr.length; i++) {
        writer.printf("%s ", inv_map.get(i));
	for (int j = 0; j < arr[i].length; j++) {
	  writer.printf("%f ", arr[i][j]);
	}
	if (curNodes.contains(i)) {
	  writer.printf("-1 ");
	}
        writer.printf("\n");
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void
  output_0d(List<Double> arr, String fileDir) {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
      for (int i = 0; i < arr.size(); i++) {
        writer.printf("%f\n", arr.get(i));
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /* arr_s.get(t) is a n*1 array */
  public static void
  output_2d(List<double[][]> arr_s, String fileDir) {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
      for (int i = 0; i < arr_s.size(); i++) {
	double[][] arr = arr_s.get(i);
        writer.printf("%d ", i);
	for (int j = 0; j < arr.length; j++) {
	  writer.printf("%f ", arr[j][0]);
	}
        writer.printf("\n");
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /* arr_s.get(t) is a n*1 array */
  public static void
  output_1d(List<double[]> arr_s, String fileDir) {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
      for (int i = 0; i < arr_s.size(); i++) {
	double[] arr = arr_s.get(i);
        writer.printf("%d ", i);
	for (int j = 0; j < arr.length; j++) {
	  writer.printf("%f ", arr[j]);
	}
        writer.printf("\n");
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /* arr_s.get(t) is a n*1 array */
  public static void
  output(Map<Integer, Integer> id_map, String fileDir) {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
      for (Map.Entry<Integer, Integer> e: id_map.entrySet()) {
	int globalID = e.getKey();
	int localID = e.getValue();
	writer.printf("%d %d\n", globalID, localID);
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }


  public static void readInit(String filename, Map<String, Integer> map, double[][] x) {
    try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
      String currentLine;
      while ((currentLine = br.readLine()) != null) {
	String[] tokens = currentLine.split(" ");
	int n = map.get(tokens[0]);
	double v1 = Double.parseDouble(tokens[1]);
	double v2 = Double.parseDouble(tokens[2]);

	x[n][0] = v1; x[n][1] = v2;
      }
    } catch (IOException _e) {
      //System.out.printf("[Info] File %s does not exist (readInit).\n", filename);
    } 
  }

}
