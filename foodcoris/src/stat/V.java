package stat;

import java.util.Iterator;
import java.util.Vector;

public class V {


    public static double proportion(double Haverage, double minH, double maxH, double minP, double maxP) {
        return ((Haverage - minH) * (maxP - minP) / (maxH - minH)) + minP;
    }

    public static Vector<Double> div(Vector<Double> vd1, Vector<Double> vd2) {
        Iterator<Double> i1 = vd1.iterator();
        Iterator<Double> i2 = vd2.iterator();
        Vector<Double> v = new Vector<Double>();
        while (i1.hasNext() & i2.hasNext()) {
            v.add(i1.next() / i2.next());
        }
        return v;
    }

    public static Vector<Double> mul(Vector<Double> vd1, Vector<Double> vd2) {
        Iterator<Double> i1 = vd1.iterator();
        Iterator<Double> i2 = vd2.iterator();
        Vector<Double> v = new Vector<Double>();
        while (i1.hasNext() & i2.hasNext()) {
            v.add(i1.next() * i2.next());
        }
        return v;
    }

    public static Vector<Double> mul(Vector<Double> vd1, Double c) {
        Iterator<Double> i1 = vd1.iterator();
        Vector<Double> v = new Vector<Double>();
        while (i1.hasNext()) {
            v.add(i1.next() * c);
        }
        return v;
    }


    public static Vector<Double> plus(Vector<Double> vd1, Double c) {
        Iterator<Double> i1 = vd1.iterator();
        Vector<Double> v = new Vector<Double>();
        while (i1.hasNext()) {
            v.add(i1.next() + c);
        }
        return v;
    }

    public static Vector<Double> negate(Vector<Double> vd1) {
        return mul(vd1, -1.0);
    }

    public static Vector<Double> f_1(Vector<Double> vd1) {
        Iterator<Double> i1 = vd1.iterator();
        Vector<Double> v = new Vector<Double>();
        while (i1.hasNext()) {
            v.add(1 / i1.next());
        }
        return v;
    }

    public static Double sum_of_elem(Vector<Double> vd) {
        Iterator<Double> i = vd.iterator();
        Double d = 0.0;
        while (i.hasNext()) {
            d += i.next();
        }
        return d;
    }

    public static Vector<Double> norm(Vector<Double> vd) {
        return mul(vd, 1 / sum_of_elem(vd));
    }

    public static Vector<Double> minus(Vector<Double> vd1, Vector<Double> vd2) {
        Iterator<Double> i1 = vd1.iterator();
        Iterator<Double> i2 = vd2.iterator();
        Vector<Double> v = new Vector<Double>();
        while (i1.hasNext() & i2.hasNext()) {
            v.add(i1.next() - i2.next());
        }
        return v;
    }

    public static Vector<Double> plus(Vector<Double> vd1, Vector<Double> vd2) {
        Iterator<Double> i1 = vd1.iterator();
        Iterator<Double> i2 = vd2.iterator();
        Vector<Double> v = new Vector<Double>();
        while (i1.hasNext() & i2.hasNext()) {
            v.add(i1.next() + i2.next());
        }
        return v;
    }

    public static Double max(Vector<Double> vd) {
        Iterator<Double> i = vd.iterator();
        Double max = Double.NEGATIVE_INFINITY;
        while (i.hasNext()) {
            max = Math.max(i.next(), max);
        }
        return max;
    }

    public static int max_num(Vector<Double> vd) {
        Iterator<Double> i = vd.iterator();
        Double max = Double.NEGATIVE_INFINITY;
        double temp = 0;
        int num = 0;
        int num_max = 0;
        while (i.hasNext()) {
            temp = i.next();
            if (temp > max) {
                max = temp;
                num_max = num;
            }
            ;
            num++;
        }
        return num_max;
    }

    public static Vector<Double> vectorDoubleIndex(Vector<Integer> vIndex, Vector<Double> vD) {
        Iterator<Integer> iI = vIndex.iterator();
        Vector<Double> fv = new Vector<Double>();
        while (iI.hasNext()) {
            fv.add(vD.elementAt(iI.next()));
        }
        return fv;
    }

    public static Vector<Double> div(double[] b, double[] c) {
        Vector<Double> vd = new Vector<Double>();
        for (int i = 0; i < Math.min(b.length, c.length); i++) {
            vd.add(b[i] / c[i]);
        }
        return vd;
    }

    public static double[] divPrime(double[] b, double[] c) {
        for (int i = 0; i < Math.min(b.length, c.length); i++) {
            if (c[i] == 0) {
                b[i] = 1; //если целевая потребность ноль, то относительная,
                //которая потом будет умножена на полное значение, занулится.
                //ставим 1, т.к. считаем, что потребность закрыта полностью
            } else {
                b[i] = b[i] / c[i];
            }
        }
        return b;
    }

    public static double[] minus(double[] b, double[] c) {
        for (int i = 0; i < Math.min(b.length, c.length); i++) {
            b[i] = b[i] - c[i];
        }
        return b;
    }

    public static double[] mul(double[] b, double[] c) {
        for (int i = 0; i < Math.min(b.length, c.length); i++) {
            b[i] = b[i] * c[i];
        }
        return b;
    }

    public static double[] plus(double[] b, double[] c) {
        for (int i = 0; i < Math.min(b.length, c.length); i++) {
            b[i] = b[i] + c[i];
        }
        return b;
    }

}
