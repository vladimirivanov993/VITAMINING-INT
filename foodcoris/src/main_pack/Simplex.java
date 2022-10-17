package main_pack;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;
import java.util.Vector;


public class Simplex {
    private static final double EPSILON = 1.0E-5;
    private double[][] a;   // tableaux
    private int M;          // number of constraints
    private int N;          // number of original variables

    private int[] basis;    // basis[i] = basic variable corresponding to row i
    // only needed to print out solution, not book

    // sets up the simplex tableaux
    public Simplex(double[][] A, double[] b, double[] c) {
        M = b.length;
        N = c.length;
        a = new double[M + 1][N + M + 1];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                a[i][j] = A[i][j];
        for (int i = 0; i < M; i++) a[i][N + i] = 1.0;
        for (int j = 0; j < N; j++) a[M][j] = c[j];
        for (int i = 0; i < M; i++) a[i][M + N] = b[i];

        basis = new int[M];
        for (int i = 0; i < M; i++) basis[i] = N + i;

        solve();

        // check optimality conditions
        assert check(A, b, c);
    }

    // run simplex algorithm starting from initial BFS
    private void solve() {
        while (true) {

            // find entering column q
            int q = bland();
            if (q == -1) break;  // optimal

            // find leaving row p
            int p = minRatioRule(q);
            if (p == -1) throw new ArithmeticException("Linear program is unbounded");

            // pivot
            pivot(p, q);

            // update basis
            basis[p] = q;
        }
    }

    // lowest index of a non-basic column with a positive cost
    private int bland() {
        for (int j = 0; j < M + N; j++)
            if (a[M][j] > 0) return j;
        return -1;  // optimal
    }

    // index of a non-basic column with most positive cost
    private int dantzig() {
        int q = 0;
        for (int j = 1; j < M + N; j++)
            if (a[M][j] > a[M][q]) q = j;

        if (a[M][q] <= 0) return -1;  // optimal
        else return q;
    }

    // find row p using min ratio rule (-1 if no such row)
    private int minRatioRule(int q) {
        int p = -1;
        for (int i = 0; i < M; i++) {
            if (a[i][q] <= 0) continue;
            else if (p == -1) p = i;
            else if ((a[i][M + N] / a[i][q]) < (a[p][M + N] / a[p][q])) p = i;
        }
        return p;
    }

    // pivot on entry (p, q) using Gauss-Jordan elimination
    private void pivot(int p, int q) {

        // everything but row p and column q
        for (int i = 0; i <= M; i++)
            for (int j = 0; j <= M + N; j++)
                if (i != p && j != q) a[i][j] -= a[p][j] * a[i][q] / a[p][q];

        // zero out column q
        for (int i = 0; i <= M; i++)
            if (i != p) a[i][q] = 0.0;

        // scale row p
        for (int j = 0; j <= M + N; j++)
            if (j != q) a[p][j] /= a[p][q];
        a[p][q] = 1.0;
    }

    // return optimal objective value
    public double value() {
        return -a[M][M + N];
    }

    // return primal solution vector
    public double[] primal() {
        double[] x = new double[N];
        for (int i = 0; i < M; i++)
            if (basis[i] < N) x[basis[i]] = a[i][M + N];
        return x;
    }

    // return dual solution vector
    public double[] dual() {
        double[] y = new double[M];
        for (int i = 0; i < M; i++)
            y[i] = -a[M][N + i];
        return y;
    }


    // is the solution primal feasible?
    private boolean isPrimalFeasible(double[][] A, double[] b) {
        double[] x = primal();

        // check that x >= 0
        for (int j = 0; j < x.length; j++) {
            if (x[j] < 0.0) {
                System.out.println("x[" + j + "] = " + x[j] + " is negative");
                return false;
            }
        }

        // check that Ax <= b
        for (int i = 0; i < M; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            if (sum > b[i] + EPSILON) {
                System.out.println("not primal feasible");
                System.out.println("b[" + i + "] = " + b[i] + ", sum = " + sum);
                return false;
            }
        }
        return true;
    }

    // is the solution dual feasible?
    private boolean isDualFeasible(double[][] A, double[] c) {
        double[] y = dual();

        // check that y >= 0
        for (int i = 0; i < y.length; i++) {
            if (y[i] < 0.0) {
                System.out.println("y[" + i + "] = " + y[i] + " is negative");
                return false;
            }
        }

        // check that yA >= c
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int i = 0; i < M; i++) {
                sum += A[i][j] * y[i];
            }
            if (sum < c[j] - EPSILON) {
                System.out.println("not dual feasible");
                System.out.println("c[" + j + "] = " + c[j] + ", sum = " + sum);
                return false;
            }
        }
        return true;
    }

    // check that optimal value = cx = yb
    private boolean isOptimal(double[] b, double[] c) {
        double[] x = primal();
        double[] y = dual();
        double value = value();

        // check that value = cx = yb
        double value1 = 0.0;
        for (int j = 0; j < x.length; j++)
            value1 += c[j] * x[j];
        double value2 = 0.0;
        for (int i = 0; i < y.length; i++)
            value2 += y[i] * b[i];
        if (Math.abs(value - value1) > EPSILON || Math.abs(value - value2) > EPSILON) {
            System.out.println("value = " + value + ", cx = " + value1 + ", yb = " + value2);
            return false;
        }

        return true;
    }

    private boolean check(double[][] A, double[] b, double[] c) {
        return isPrimalFeasible(A, b) && isDualFeasible(A, c) && isOptimal(b, c);
    }

    // print tableaux
    public void show() {
        System.out.println("M = " + M);
        System.out.println("N = " + N);
        for (int i = 0; i <= M; i++) {
            for (int j = 0; j <= M + N; j++) {
                System.out.printf("%7.2f ", a[i][j]);
            }
            System.out.println();
        }
        System.out.println("value = " + value());
        for (int i = 0; i < M; i++)
            if (basis[i] < N) System.out.println("x_" + basis[i] + " = " + a[i][M + N]);
        System.out.println();
    }


    public static void test(double[][] A, double[] b, double[] c) {
        Simplex lp = new Simplex(A, b, c);
        System.out.println("value = " + lp.value());
        double[] x = lp.primal();
        for (int i = 0; i < x.length; i++)
            System.out.println("x[" + i + "] = " + x[i]);
        double[] y = lp.dual();
        for (int j = 0; j < y.length; j++)
            System.out.println("y[" + j + "] = " + y[j]);
    }

    public static void test1() {
        double[][] A = {
                {-1, 1, 0},
                {1, 4, 0},
                {2, 1, 0},
                {3, -4, 0},
                {0, 0, 1},
        };
        double[] c = {1, 1, 1};
        double[] b = {5, 45, 27, 24, 4};
        test(A, b, c);
    }


    // x0 = 12, x1 = 28, opt = 800
    public static void test2() {
        double[] c = {13.0, 23.0};
        double[] b = {480.0, 160.0, 1190.0};
        double[][] A = {
                {5.0, 15.0},
                {4.0, 4.0},
                {35.0, 20.0},
        };
        test(A, b, c);
    }

    // unbounded
    public static void test3() {
        double[] c = {2.0, 3.0, -1.0, -12.0};
        double[] b = {3.0, 2.0};
        double[][] A = {
                {-2.0, -9.0, 1.0, 9.0},
                {1.0, 1.0, -1.0, -2.0},
        };
        test(A, b, c);
    }

    // degenerate - cycles if you choose most positive objective function coefficient
    public static void test4() {
        double[] c = {10.0, -57.0, -9.0, -24.0};
        double[] b = {0.0, 0.0, 1.0};
        double[][] A = {
                {0.5, -5.5, -2.5, 9.0},
                {0.5, -1.5, -0.5, 1.0},
                {1.0, 0.0, 0.0, 0.0},
        };
        test(A, b, c);
    }


    public static void test5() {
        double[] c = {32, 8};
        double[] b = {1, 7};
        double[][] A = {
                {1, -1},
                {3, 1},
        };
        test(A, b, c);
    }

    public static void test6() {
        double[] c = {12, 9};
        double[] b = {63, 147, 126};
        double[][] A = {
                {9, 3},
                {7, 21},
                {9, 10},
        };
        test(A, b, c);
    }

    public static void test7() {
        double[] b = {4, 3};
        double[] c = {63, 217, 126};
        double[][] A = {
                {9, 7, 9},
                {3, 21, 10}
        };
        test(A, b, c);
    }

    public static void test8() {
        double[] b = {0.04, 0.15, 0.40};
        double[] c = {0.8, 22, -5};
        double[][] A = {
                {0.38, 0, 0},
                {0.001, 0.09, -0.02},
                {0.002, 0.5, -0.8}
        };
        test(A, b, c);
    }

    public static void test9() {
        double[] b = {450, 400};
        double[] c = {45, 80};
        double[][] A = {
                {10, 15},
                {5, 20},
        };
        test(A, b, c);
    }


    public static void test10() {
        double[] b = {1, 1, 1, 1, 1, 1, 1};//,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1, 1, 1,1,   1,1,1,1,  1,1,1,1,1};
        //double[] c = {75, 0, 1.5, 2.0, 20, 0.25, 20, 1500, 2000, 1125, 400, 1000, 19, 0.1};
        double[] c = {75, 0, 1.5, 2.0, 20, 0.0, 0.0, 1500, 2000, 1125, 400, 1000, 19, 0.0};
        System.out.println("b=" + b.length);


        double[][] A = {{70, 0.0, 0.02, 0.1, 0.6, 0.0, 0.0, 10, 210, 26, 17, 51, 1, 0.0
        }, {10, 0.0, 0.06, 0.04, 0.2, 0.0, 0.0, 8, 141, 23, 14, 42, 1, 0.0
        }, {25, 0.0, 0.0, 0.04, 0.1, 0.0, 0.0, 10, 255, 39, 13, 44, 1, 0.0
        }, {15, 0.0, 1.75, 0.08, 0.05, 0.0, 0.0, 8, 220, 77, 40, 34, 1, 0.0
        }, {25, 0.0, 1.2, 0.04, 0.53, 0.0, 0.0, 40, 290, 14, 20, 26, 1, 0.0
        }, {15, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 7, 73, 40, 7, 16, 0.4, 0.0
        }, {2.8, 0.0, 0.05, 0.1, 1, 0.0, 0.0, 3, 664, 124, 198, 564, 2, 0.0
        }};

        System.out.println("A=" + A.length);

        test(A, b, c);
    }


    public static double[][] M1() throws FileNotFoundException {
        Vector<String> vs1 = stat.M.translatorFromFile(new File("in.txt"));
        return stat.M.extraxtMatrix(vs1, 1);
    }

    ;


    public static void test11() throws FileNotFoundException {
//      double[] b = {1,1, 1,1,   1,1,1      };
        double[] b = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        System.out.println("b=" + b.length);
        double[] c = {75, 0, 1.5, 1.5, 2.0, 20, 2, 0.05, 0.4, 0.25, 0.2, 15, 20, 1500, 2000, 1125, 400, 1000, 19, 0.1};
        System.out.println("c=" + c.length);
        //double[][] A = M1();
        Vector<String> vs = stat.M.translatorFromFile(new File("in.txt"));
        double[][] A = stat.M.extraxtMatrix(vs, 1);
        Vector<Integer> vi1 = new Vector<Integer>();
        vi1.addElement(0);
        Vector<Integer> vi2 = new Vector<Integer>();
        vi2.addElement(0);
        //A = stat.M.extraxtPartMatrixModl(A, vi1, vi2);
        System.out.println("A n=" + A.length);
        double[] am = A[0];
        double[] b1 = {1};
        double[] c1 = {75};
        System.out.println("A m=" + am.length);
        test(A, b, c);

        Simplex lp = new Simplex(A, b, c);//lp.primal;lp.primal;lp.value();
        Vector<String> names = stat.M.positiveName(vs, lp.dual());
        stat.M.printToFile(names, "names1.txt");
        //Vector<String> names2 = stat.M.positiveNameWithMean(vs, lp.dual());
        //stat.M.printToFile(names2, "names2.txt");


        stat.M.FileConvertBZHU(new File("продукты бжу.txt"), new File("bzhu.txt"));
        Vector<String> vs3 = stat.M.translatorFromFile(new File("bzhu.txt"));
        double[][] A3 = stat.M.extraxtMatrix(vs3, 1);
        double[] b3 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        double[] c3 = {500, 162, 0, 0, 00};
        Simplex lp3 = new Simplex(A3, b3, c3);
        System.out.println("begin");
        test(A3, b3, c3);
        System.out.println("end");
        System.out.println("simplex_val=" + lp3.value());
        Vector<String> names3 = stat.M.positiveName(vs3, lp3.dual());
        stat.M.printToFile(names3, "names3.txt");


    }


    public static void test12() throws FileNotFoundException {
//    double[] b = {1,1, 1,1,   1,1,1      };
        double[] b = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        System.out.println("b=" + b.length);
//		double[] c = {500, 162,113.4,0, 2400,75,0,1.5,1.5,2.0,20,2,0.05,0.4,0.25,0.2,15,20,1500,2000,1125,400,1000,19,0.1};
        double[] c = {500, 120, 72, 480, 2400, 75, 0, 1.5, 1.5, 2.0, 20, 2, 0, 0.4, 0.25, 0.2, 0, 20, 1500, 2000, 1125, 400, 1000, 19, 0.1};
        System.out.println("c=" + c.length);
        Vector<String> vs = stat.M.translatorFromFile(new File("in5.txt"));
        double[][] A = stat.M.extraxtMatrix(vs, 1);
        Simplex lp = new Simplex(A, b, c);//lp.primal;lp.primal;lp.value();
        Vector<String> names = stat.M.positiveName(vs, lp.dual());
        stat.M.printToFile(names, "names5.txt");


    }


    public static void test13() throws FileNotFoundException {
//    double[] b = {1,1, 1,1,   1,1,1      };
        double[] b = {23, 100, 72, 439, 200, 19, 52, 20, 15, 17, 24, 250, 40, 35, 75, 235, 95, 110, 520, 537, 50, 460, 115, 450, 162, 60, 75, 60, 48, 180, 65, 79, 440, 45, 65, 75, 103, 350, 260, 168, 90, 165, 200, 51, 98, 85, 180, 54, 65, 59, 75, 40, 170, 173, 119, 400, 43, 26};
        System.out.println("b=" + b.length);
//		double[] c = {500, 162,113.4,0, 2400,75,0,1.5,1.5,2.0,20,2,0.05,0.4,0.25,0.2,15,20,1500,2000,1125,400,1000,19,0.1};
        double[] c = {100, 80, 65, 320, -6000, 75, 0, 1.5, 1.5, 2.0, 20, 2, 0, 0.4, 0.25, 0.2, 0, 20, 1500, 2000, 1125, 400, 1000, 19, 0.1};
        System.out.println("c=" + c.length);
        Vector<String> vs = stat.M.translatorFromFile(new File("in5.txt"));
        double[][] A = stat.M.extraxtMatrix(vs, 1);
        Simplex lp = new Simplex(A, b, c);//lp.primal;lp.primal;lp.value();
        Vector<String> vs2 = new Vector<String>();
        vs2.add("вода");
        vs2.add("белки");
        vs2.add("жиры");
        vs2.add("углеводы");
        vs2.add("калории");
        vs2.add("vitamin_C");
        vs2.add("vitamin_P");
        vs2.add("vitamin_A");
        vs2.add("vitamin_B1");
        vs2.add("vitamin_B2");
        vs2.add("vitamin_PP");
        vs2.add("vitamin_B6");
        vs2.add("vitamin_B12");
        vs2.add("vitamin_Bc");
        vs2.add("vitamin_K");
        vs2.add("vitamin_H");
        vs2.add("vitamin_D");
        vs2.add("vitamin_E");
        vs2.add("mineral_Na");
        vs2.add("mineral_K");
        vs2.add("mineral_Ca");
        vs2.add("mineral_Mg");
        vs2.add("mineral_P");
        vs2.add("mineral_Fe");
        vs2.add("mineral_J");
        double[] c_constr = {150, 10, 100, 100, 30, 150, 200, 200, 150, 150, 25, 100, 100, 100, 100, 70, 500, 500, 20, 70, 70, 30, 30, 30, 30, 150, 150, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 15, 5, 200, 180};
        c_constr = stat.M.mulC(c_constr, 0.01);
        double[] c_realNext = c_constr;

        c_realNext = i_entry_simplex(A, b, c, lp, vs, vs2, c_constr, "names5_1.txt");//возвращаем проценты от нормы витамин
        double[] y_int = c_realNext; //y1 - добавляем данные, полученные за первый день
        double n = 2.0; //считаем для второго дня - должны были набрать 200%, а в абсолюном значении это удвоенная норма
        int n_int = (int) n;
        double[] c_realNextCorrect = stat.M.minus_overload_iterable(y_int, n); //если набрали в переизбытке, учитываем его. на выходе процент от нормы по витаминам, который необходимо будет набрать на следующем шаге
        c_realNext = stat.V.mul(c_realNextCorrect, c); //из относительного значения переводим в абсолютное

        while (n < 3) { //указываем, на сколько дней отправить вперёд
            lp = new Simplex(A, b, c_realNext);//lp.primal;lp.primal;lp.value();
            c_realNext = i_entry_simplex(A, b, c_realNext, lp, vs, vs2, c_constr, "names5_" + n_int + ".txt");//возвращаем проценты от нормы витамин
            y_int = stat.V.plus(y_int, c_realNext); //y1 - добавляем данные, полученные за первый день
            n += 1.0; //считаем для второго дня - должны были набрать 200%, а в абсолюном значении это удвоенная норма
            n_int = (int) n;
            c_realNextCorrect = stat.M.minus_overload_iterable(y_int, n); //если набрали в переизбытке, учитываем его. на выходе процент от нормы по витаминам, который необходимо будет набрать на следующем шаге
            c_realNext = stat.V.mul(c_realNextCorrect, c); //процент по витаминам переводим из относительного значения в абсолютное. и данные для нового симплекса готовы
        }


        //матрица вкладов по витаминам(умножаем каждый слолбец оригинальной матрицы на количество продукты)
        //новый день: вычитаем дневную потребность из переизбытка(переизбыток не может сменить знак своего ограничения)
        //и по тем витаминам, по которым ещё не закрыли потребность, закрываем
        //алгоритм мчится на 1 месяц вперед
        //после того, как формировали потребность по витаминам на месяц, планируем закупки с учетом срока годности
        //продумать момент, что если потребность по витамину превышена в 10 раз, то ограничения переизбытка не более 10ти
        //-идея относительно отрицательных значений - они не учитываются, их переизбыток всегда ноль
        //вообще, перед оптимизацией витамин, нормировать их, чтобы сразу получать результат в процентах и не путаться
        //потом это пригодится для оптимизации
        //некоторые витамины плохо накапливаются - подумать, стоит ли их включать в вектор переизбытка
        //переизбыток по итогам оптимизации может быть любым. но после 1 месяца(или заданного интервала) считаем его нулевым
        //виновник переизбытка - продукт, общее количество которого в показателе по витамину многократно превышает его норму

        //в оптимизации можно выбирать подмножества продуктов, и запускать только для них.
        //тут подойдет нейросеть - количество выходных нейронов может быть по количеству ограничивающих продутов
        //прогамма оптимизации будет давать числовой результат для каждого набора продуктов
        //по сути, чем больше множество продуктов, тем точнее будет выбор.
        //но если мы заранее не хотим что-то, можем это убрать
        //нейросеть можно обучить вручную - понравилось - не понравилось(множество тех продуктов, которые попали в результат)
        //по результатам обучения нейронка сама сможет определять, понравится ли будующий набор
        //у продуктов будет статистика, как часто он понравился или не понравился, это будет влиять на сужение оптимизационного множества
        //программа будет стартовать с разных начальных условий, и то, которое оценит, как потенциально понравившиеся, выдаст
        //при первом запуске будет предложено выбрать, какие продукты нравятся больше всего(выбрать множество оптимизации)
    }


    //полученные разложения нужно сложить

    /**/


    public static double[] i_entry_simplex(double[][] A, double[] b, double[] c, Simplex lp,
                                           Vector<String> vs, Vector<String> vs2, double[] c_constr, String filename) throws FileNotFoundException {
        Vector<String> names = stat.M.positiveName(vs, lp.dual()); //имена продуктов-победителей и их количествах, пойдут в общий лог 
        double[] c_real = stat.M.matrixTransposeMulVect(A, lp.dual());//количество витамин, собраных в продуктах: общее количество каждого витамина по всем продуктам
        Vector<String> names2 = stat.M.positiveName(vs2, c_real);  //количества витамин, в результирующей выборке
        names.addAll(names2);                                       //добавим в лог вывода собранные витамины
        names.add("cost=" + Double.toString(lp.value()));               //выведем стоимость набора продуктов в лог
        double[] ar = stat.M.razlozheniePoOgranicheniu(lp.dual(), c_constr)[0];//оставшиеся количество продуктов, после вычета излишков
        double[] lr = stat.M.razlozheniePoOgranicheniu(lp.dual(), c_constr)[1];//все излишки по продуктам
        Vector<Integer> vinn = stat.M.notZero(lr);                             //выясним, по каким номерам продуктов нашли излишки
        double[] res;
        if (vinn.size() != 0) {
            double[][] At = stat.M.matrixTranspose(A);   //транспонируем матрицу
            double[] cNew = stat.M.matrixMulVect(At, ar);//получаем вектор из витамин по продуктам, оставшихся после вычета излишков, элемент вектора = сумма витамина по всем продуктам
            double[][] A0 = stat.M.zeroid(A, vinn);      //мы зануляем количества витамин по тем продуктам, по которым были излишки, чтобы убрать их из оптимизации
            double[] cOst = stat.M.minusCPost(c, cNew);  //считаем, сколько витамин нам теперь недостаёт
            Simplex lpNew;
            double[] lpNextY;
            try {
                lpNew = new Simplex(A0, b, cOst);     //добираем недостающие;
                lpNextY = stat.M.plus(ar, lpNew.dual()); //добавляем к нашему набору новые продукты, чтобы восполнить потерянные витамины
            } catch (java.lang.ArithmeticException e) {
                System.out.println("разложение по ограничению не работает - симлекс не имеет решения");
                lpNew = lp;
                lpNextY = ar; //сколько уже получилось, больше не собрать, не нарушив условия
            }
            double[] c_realNext = stat.M.matrixMulVect(At, lpNextY); //получаем вектор полученных витамин в продуктовом наборе
            Vector<String> names3 = stat.M.positiveName(vs, lpNextY);//выясняем имена этих продуктов
            names.addAll(names3);                                    //отправляем в лог печати
            Vector<String> names4 = stat.M.positiveNameProcent(vs2, c_realNext, c);//получаем абсолютное кол-во витамин и %
            names.addAll(names4);                                     //отправляем в лог печати
            names.add("cost=" + Double.toString(lpNew.value()));         //считаем стоимость нового набора
            res = c_realNext; //возвращаем количество полученых витамин
        } else {
            res = c_real; //возвращаем количество полученых витамин
        }
        res = stat.V.divPrime(res, c); //получаем процентное отношение
        stat.M.printToFile(names, filename);
        return res;
    }
	/*
    
    public static void test5() {
        double[] c = { 32, 8};
        double[] b = {  1, 7 };
        double[][] A = {
            { 1, -1},
            { 3, 1},
        };
        test(A, b, c);
    }
    
    */

    public static void test0000(double[][] A, double[] b, double[] c) {
        Simplex lp = new Simplex(A, b, c);
        System.out.println("value = " + lp.value());
        double[] x = lp.primal();
        for (int i = 0; i < x.length; i++)
            System.out.println("x[" + i + "] = " + x[i]);
        double[] y = lp.dual();
        for (int j = 0; j < y.length; j++)
            System.out.println("y[" + j + "] = " + y[j]);
    }


    // test client
    public static void main(String[] args) throws FileNotFoundException {

        /*try                           { test1();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                           { test2();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                           { test3();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                           { test4();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                           { test5();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                           { test6();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                           { test7();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");        

        try                           { test8();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");
        
        try                           { test9();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");
        
        
        try                           { test10();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");   
        
        
        try                           { test11();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");
        
        try                           { test12();             }
        catch (ArithmeticException e) { e.printStackTrace(); }
        System.out.println("--------------------------------");*/

        try {
            test13();
        } catch (ArithmeticException e) {
            e.printStackTrace();
        }
        System.out.println("--------------------------------");

        
        /*
        int M = Integer.parseInt(args[0]);
        int N = Integer.parseInt(args[1]);
        double[] c = new double[N];
        double[] b = new double[M];
        double[][] A = new double[M][N];
        for (int j = 0; j < N; j++)
            c[j] = StdRandom.uniform(1000);
        for (int i = 0; i < M; i++)
            b[i] = StdRandom.uniform(1000);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A[i][j] = StdRandom.uniform(100);
        Simplex lp = new Simplex(A, b, c);
        System.out.println(lp.value());
        */
    }
}



