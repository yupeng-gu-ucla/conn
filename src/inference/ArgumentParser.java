import java.io.*;
import java.util.*;

/* options:
 *   -data: the directory of data files
 *   -sigma: value of sigma; 1.0 by default
 *   -eps: value of eps; 0.8 by default
 *   -iter: value of maximum iterations for stochastic gradient ascent (in thousands); 5000 thousand by default
 *   -outer: value of maximum outer iterations for coordinate gradient ascent; 7 by default
 *   -verbose: whether see verbose output or not (0: show all outputs; other: show limited outputs); 1 by default
 *   -start: start timeslice; 120 by default
 *   -time: length of timespan; 100 by default
 */

public class ArgumentParser {
  public static void parse(String[] args) {
    int k = 0, K = args.length;
    while (k < K) {
      switch (args[k]) {
	case "-sigma":
	  Main.sigma = Double.parseDouble(args[++k]);
	  break;
	case "-eps":
	  Main.eps = Double.parseDouble(args[++k]);
	  break;
	case "-iter":
	  Main.MAX_ITER = 1000 * Integer.parseInt(args[++k]);
	  break;
	case "-outer":
	  Main.NUM_OUTER_ITER = Integer.parseInt(args[++k]);
	  break;
	case "-verbose":
	  int v = Integer.parseInt(args[++k]);
	  Main.VERBOSE = (v != 0);
	  break;
	case "-data":
	  Main.prefix = args[++k];
	  break;
	case "-start":
	  Main.tStart = Integer.parseInt(args[++k]);
	  break;
	case "-time":
	  Main.T = Integer.parseInt(args[++k]);
	  break;
      }
      k++;
    }
  }

  public static void help() {
    String h = "\nUsage: java Main <options>\n" 
      + "options:\n" 
      + "\t-data: the directory of data files\n"
      + "\t-sigma: value of sigma (1.0 by default)\n"
      + "\t-eps: value of eps (0.8 by default)\n"
      + "\t-iter: value of maximum iterations for stochastic gradient ascent, in thousands (5000 by default)\n"
      + "\t-outer: value of maximum outer iterations for coordinate gradient ascent (7 by default)\n"
      + "\t-verbose: whether see verbose output or not; 0: show all outputs; other: show limited outputs (1 by default)\n"
      + "\t-start: start timeslice (120 by default)\n"
      + "\t-time: length of timespan (100 by default)\n";
    System.out.println(h);
  }
}

