public class Util {

  public static double logis(double x) {
    if (x < -20) {
      return Double.MIN_VALUE;
    } else {
      double res = 1.0 / (1.0 + Math.exp(-x));
      return (res > 1.0-1e-9) ? 1.0-1e-9 : res;
    }
  }


  /* Distance between two directions */
  /* return value between -Pi and Pi */
  public static double angleDistance(double a1, double a2) {
    while (a1 < -Math.PI) a1 += 2*Math.PI;
    while (a1 > Math.PI)  a1 -= 2*Math.PI;
    while (a2 < -Math.PI) a2 += 2*Math.PI;
    while (a2 > Math.PI)  a2 -= 2*Math.PI;
    if (Math.abs(a1-a2) < Math.PI) return a1-a2;
    else return 2*Math.PI - Math.abs(a1-a2);
  }

}


