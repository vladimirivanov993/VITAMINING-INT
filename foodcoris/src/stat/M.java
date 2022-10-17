package stat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Vector;

public class M {
    //replace space 0
    //replace \t ,
    //first name into product name
    //количество столбцов берется из первой строки
    //количество строк таблицы из количества всех строк
    //(File file1, int m, int n, int offset)
    public M() {
    }

    /**
     * matrix M x N
     *
     * @offset количество пропускаемых строк с начала. если только название, значение в 1. если больше 1, название обяз.
     */

    public static Vector<String> translatorFromFile(File file1) throws FileNotFoundException {
        Scanner sc = new Scanner(file1);
        Vector<String> matrix = new Vector<String>();
        while (sc.hasNextLine()) {
            matrix.add(sc.nextLine());//автоматическое преобразование к нужному, осавить в виде строки
        }
        sc.close();
        return matrix;
    }

    /**
     * НЕЗАВИСИМ. имя всегда должно быть первым, после него запятая
     */
    public Vector<String> getNameProduct(Vector<String> matrix) {
        Iterator<String> i = matrix.iterator();
        Vector<String> names = new Vector<String>();
        while (i.hasNext()) {
            names.addElement(i.next().split(",")[0]);
        }
        return names;
    }

    public static String convertStringForParse(String s) {
        System.out.println("1=" + s + "=!");
        //s=s.replaceAll("\r", "0.0");
        //s=s.replaceAll("\t", ";0");
        //System.out.println("2="+s+"=!");
        //s=s.replaceAll(";", ";0");
        //System.out.println("3="+s+"=!");

        return s;
    }

    public static String[] splitMeans(String s) {
        return s.split(";");
    }

    /**
     * Гречка,2/3 ст,0,0,0,0.43,0.2   offset=2
     */
    public static double[][] extraxtMatrix(Vector<String> vs, int offset) {

        String s1 = convertStringForParse(vs.firstElement());//заголовок таблицы и названия витаминов
        String[] ms1 = s1.split(";");
        int M = Math.max(0, ms1.length - offset);
        int N = Math.max(0, vs.size());
        double[][] matrix = new double[N][M];
        Iterator<String> i = vs.iterator();
        String[] temp = new String[M];
        int indexN = 0;
        while (i.hasNext()) {
            temp = splitMeans(convertStringForParse(i.next()));
            for (int k = 0; k < M; k++) {
                matrix[indexN][k] = Double.parseDouble(temp[k + offset]);
            }
            indexN++;
        }
        return matrix;
    }


    /**
     * vM и vN должны содержать множество числовых значений. порядок создания матрицы определяется ими же
     * НЕЗАВИСИМ. другие определяют лишь правильную структуру
     */
    public static Double[][] extraxtPartMatrix(Double[][] matrix, Vector<Integer> vM, Vector<Integer> vN) {
        Double[][] retMatrix = new Double[vN.size()][vM.size()];
        Iterator<Integer> iN = vN.iterator();
        int iNc;
        int rM = 0;
        int rN = 0;
        while (iN.hasNext()) {
            Iterator<Integer> iM = vM.iterator();
            iNc = iN.next();
            rN++;
            rM = 0;
            while (iM.hasNext()) {
                retMatrix[rN][rM] = matrix[iNc][iM.next()];
                rM++;
            }
        }
        return retMatrix;
    }

    /*выбор из матрицы множества подстрок*/
    public static double[][] extraxtPartMatrixSubRow(double[][] a, Vector<Integer> vN) {
        double[][] retMatrix = new double[vN.size()][a[0].length];
        Iterator<Integer> iN = vN.iterator();
        int rN = 0;
        while (iN.hasNext()) {
            rN++;
            retMatrix[rN] = a[iN.next()];
        }
        return retMatrix;
    }


    public static double[] subSet(double[] a, Vector<Integer> vi) {
        double[] ret = new double[vi.size()];
        Iterator<Integer> i = vi.iterator();
        int index = 0;
        while (i.hasNext()) {
            ret[index] = a[i.next()];
            index++;
        }
        return ret;
    }


    public static double[][] extraxtPartMatrixSubCol(double[][] matrix, Vector<Integer> vM) {
        double[][] retMatrix = new double[matrix.length][vM.size()];
        for (int i = 0; i < matrix.length; i++)
            retMatrix[i] = subSet(matrix[i], vM);
        return retMatrix;
    }


    /**
     * альтернативная редакция
     */
    public static double[][] extraxtPartMatrixModl(double[][] a, Vector<Integer> vM, Vector<Integer> vN) {
        double[][] retMatrix = new double[vN.size()][vM.size()];
        Iterator<Integer> iN = vN.iterator();
        int index = 0;
        while (iN.hasNext()) {
            retMatrix[index] = subSet(a[iN.next()], vM);
            index++;
        }
        return retMatrix;
    }

    /**
     * поиск ненулевых номеров
     */
    public static Vector<String> selectStringOnNotNul(Vector<String> vs, Vector<Double> vd) {
        Iterator<String> is = vs.iterator();
        int index = 0;
        Vector<String> ret = new Vector<String>();
        while (is.hasNext()) {
            if (Math.abs(vd.elementAt(index)) == 0) {
                ret.addElement(is.next());
            }
            index++;
        }
        return ret;
    }

    /**
     * выборка из массива строк по массиву целых
     */
    public static Vector<String> selectStringsOnNumbers(Vector<String> vs, Vector<Integer> vi) {
        Iterator<Integer> ii = vi.iterator();
        Vector<String> ret = new Vector<String>();
        while (ii.hasNext()) {
            ret.addElement(vs.elementAt(ii.next()));
        }
        return ret;
    }


    /**
     * поиск ненулевых строк и возвращение имени перед первым разделителем точка с запятой
     */
    public static Vector<String> positiveName(Vector<String> vs, double[] vd) {
        Vector<String> ret = new Vector<String>();
        Double d = new Double(0);
        for (int index = 0; index < Math.min(vd.length, vs.size()); index++) {
            if (Math.abs(vd[index]) > 0) {
                d = vd[index];
                ret.addElement((d.toString().substring(0, Math.min(d.toString().length(), 7))).concat("\t" + vs.elementAt(index).split(";")[0]));
            }
        }
        return ret;
    }

    /**
     * поиск ненулевых строк и возвращение имени c значением перед первым разделителем точка с запятой
     */
    public static Vector<String> positiveNameWithMean(Vector<String> vs, double[] vd) {
        Vector<String> ret = new Vector<String>();
        Double d = new Double(0);
        for (int index = 0; index < Math.min(vd.length, vs.size()); index++) {
            if (Math.abs(vd[index]) > 0) {
                d = vd[index];
                ret.addElement(vs.elementAt(index).split(";")[0].concat(d.toString()));
            }
        }
        return ret;
    }


    /**
     * получение номеров нулевых значений
     */
    public static Vector<Integer> zeroNums(double[] vd) {
        Vector<Integer> ret = new Vector<Integer>();
        for (int i = 0; i < vd.length; i++) {
            if (Math.abs(vd[i]) == 0) {
                ret.addElement(i);
            }
        }
        return ret;
    }

    /**
     * получение номеров ненулевых значений
     */
    public static Vector<Integer> notZero(double[] vd) {
        Vector<Integer> ret = new Vector<Integer>();
        for (int i = 0; i < vd.length; i++) {
            if (Math.abs(vd[i]) != 0) {
                ret.addElement(i);
            }
        }
        return ret;
    }


    public static void printVectorString(Vector<String> vs) {
        Iterator<String> is = vs.iterator();
        while (is.hasNext()) {
            System.out.println(is.next());
        }
    }

    /**
     * печать в файл по заданному имени
     */
    public static void printToFile(Vector<String> vs, String filename) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new File(filename));
        Iterator<String> is = vs.iterator();
        while (is.hasNext()) {
            pw.write(is.next() + "\n");
        }
        pw.close();
    }

    public static double[] plus(double[] a, double[] b) {
        for (int i = 0; i < Math.min(a.length, b.length); i++) {
            a[i] += b[i];
        }
        return a;
    }

    /**
     * Ax=b
     */
    public static double[] matrixMulVect(double[][] A, double[] x) {
        int M = A.length;
        int N = A[0].length;
        System.out.println("M=" + M);
        System.out.println("N=" + N);
        System.out.println("x.length=" + x.length);
        //строка на столбец
        double[] b = new double[M];
        if (A[0].length == x.length) {
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    b[i] += A[i][j] * x[j];
                }
            }
        } else return null;
        return b;
    }

    /**
     * matrix transpose
     */
    public static double[][] matrixTranspose(double[][] A) {
        double[][] At = new double[A[0].length][A.length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                //System.out.println(
                At[j][i] = A[i][j];//);
            }
        }
        return At;
    }


    public static double[] matrixTransposeMulVect(double[][] A, double[] x) {
        int M = A.length;
        int N = A[0].length;
        System.out.println("M=" + M);
        System.out.println("N=" + N);
        System.out.println("x.length=" + x.length);
        //строка на столбец
        double[] b = new double[N];
        if (M == x.length) {
            for (int j = 0; j < M; j++)
                for (int i = 0; i < N; i++) {
                    {
                        b[i] += A[j][i] * x[j];
                    }
                }
        } else return null;
        return b;
    }

    public static double[] mulC(double[] a, double c) {
        int l = a.length;
        double[] b = new double[l];
        for (int i = 0; i < l; i++) {
            b[i] = a[i] * c;
        }
        return b;
    }

    public static double[] minusC(double[] a, double[] b) {
        if (a.length == b.length) {
            int l = a.length;
            //double[] c = new double[l];
            for (int i = 0; i < l; i++) {
                b[i] = a[i] - b[i];
            }
            return b;
        } else
            return null;
    }

    public static double[] minusCPost(double[] a, double[] b) {
        if (a.length == b.length) {
            int l = a.length;
            double[] c = new double[l];
            for (int i = 0; i < l; i++) {
                c[i] = a[i] - b[i];
                if ((a[i] >= 0 || b[i] >= 0) && c[i] < 0) c[i] = 0;
            }
            return c;
        } else
            return null;
    }

    /**
     * @l - вектор ограничения
     * @a - вектор исходного распределения
     */
    public static double[][] razlozheniePoOgranicheniu(double[] a, double[] l) {
        double[][] razl = new double[2][a.length];
        if (a.length == l.length) {
            for (int i = 0; i < a.length; i++) {
                //if(Math.abs(a[i]-l[i]) > 0.05*Math.abs(l[i])) 14янв
                if ((a[i] - l[i]) > 0.05 * l[i]) {
                    razl[0][i] = l[i];
                    razl[1][i] = a[i] - l[i];
                } else {
                    razl[0][i] = a[i];
                    razl[1][i] = 0.0;
                }

            }
        } else return null;
        return razl;
    }


    /**
     * сложение одгого вектора с другим по индексу второго из элементов первого
     * индекс второго не должен превосходить количества элементов первого
     */
    public static double[] plusIndex(double[] a, double[] b, Vector<Integer> indexB) {
        double[] ret = a;
        Iterator<Integer> vi = indexB.iterator();
        int i = 0;
        while (vi.hasNext()) {
            ret[vi.next()] += b[i];
            i++;

        }
        return ret;


    }


    public static Vector<Integer> getPositiveNumbers(double[] a) {
        Vector<Integer> vi = new Vector<Integer>();
        for (int i = 0; i < a.length; i++) {
            if (a[i] > 0) vi.add(i);
        }
        return vi;
    }


    //ПОКА НЕ ИСПОЛЬЗУЕТСЯустранение из задачи оптимизации поднабора витамин,
    //ограничения по которым заведомо не могут быть выполнены
    //строка содержит распределение витаминов по продукту
    //на выходе вектор содержит порядковые номера витамин

    /***/
    public static double[] ixcludeNullVitamin(double[][] A) {
        double[] a = new double[A.length];
        for (int i = 0; i < A.length; i++) {
            a = plus(a, A[0]);
        }
        return a;
    }


    /***/
    //выбор меры с ненулевым содержанием по витаминам.
    //по каждому витамину определяется мера, определяющая, в каком количестве продуктов из перечисленных он содержится
    public static Vector<Double> extractMeasuresVita() {
        return null;

    }

    public static void FileConvertBZHU(File in, File out) throws FileNotFoundException {
        Scanner sc = new Scanner(in);
        String s = "";
        String names = "";
        boolean header = false;
        String prefix = "";
        PrintWriter pw = new PrintWriter(out);
        while (sc.hasNextLine()) {
            //System.out.println("iamlive!");
            s = sc.nextLine();
            if (!(s.contains(";"))) {
                //header = true;
                if (!s.equals("")) prefix = s;
                System.out.println("prefix=" + prefix);
            }
            if (s.contains(";")) {
                //s.replaceAll(",", ".");
                pw.write(prefix.concat("-").concat(s) + "\n");
                System.out.println(prefix.concat("-").concat(s) + "\n");
                header = false;
            }
        }
        sc.close();
        pw.close();
    }

    /*зануление заданных строк в матрице*/
    public static double[][] zeroid(double[][] A, Vector<Integer> vinn) {
        // TODO Auto-generated method stub
        Iterator<Integer> vi = vinn.iterator();
        int ind = 0;
        while (vi.hasNext()) {
            ind = vi.next();
            for (int i = 0; i < A[0].length; i++) {
                A[ind][i] = 0;
            }
        }
        return A;
    }

    public static Vector<String> positiveNameProcent(Vector<String> vs2, double[] b_realNext, double[] b) {
        // TODO Auto-generated method stub
        Vector<String> pn = positiveName(vs2, b_realNext);
        Iterator<String> is = pn.iterator();
        //74.9999	C
        String s;
        String s_pn;
        int vitaindex = 0;
        Double d_percent;
        Vector<String> vs_res = new Vector<String>();
        while (is.hasNext()) {
            s_pn = is.next();
            s = s_pn.split("\t")[1]; //имя витамина
            vitaindex = vs2.indexOf(s);
            d_percent = (b_realNext[vitaindex] / b[vitaindex]) * 100;
            s = s_pn + "\t(" + d_percent.toString().substring(0, Math.min(d_percent.toString().length(), 7)) + "%)";
            vs_res.add(s);
        }
        return vs_res;

    }

    /**
     * если какое-то значение у вектора отрицательное, его не меняем.
     * предполагаем, что на вход подаются доли от целого(может и больше)
     * а на выходе определяем, сколько нужно добрать
     */
    public static double[] minus_overload(double[] c) {
        // TODO Auto-generated method stub
        for (int i = 0; i < c.length; i++) {
            if (c[i] > 0.99) {//если прошлый раз не добрали немного - это погрешность, т.к. по цели задачи должны были покрыть в пределах погрешности
                // c = 1.4 - было превышение на 40%, должно быть 0.6 в итоге
                // с = 1.0 - прошлый раз выработали в 100%, с не меняется, цель снова 100%
                // с < 0.99.. - немного не добрали в прошлый раз, но цель была добрать максимум, поэтому не меняем ничего - это погрешность
                // c > 2.1 - запасов хватит на три дня или более, нет никакой потребности
                c[i] = Math.max(0, 2.0 - c[i]);
            }
        }
        return c;
    }

    public static double[] minus_overload_iterable(double[] y_int, double n) {
        for (int i = 0; i < y_int.length; i++) {
            if (y_int[i] > 0) {//если прошлый раз не добрали немного - это погрешность, т.к. по цели задачи должны были покрыть в пределах погрешности
                // c = 1.4 - было превышение на 40%, должно быть 0.6 в итоге
                // с = 1.0 - прошлый раз выработали в 100%, с не меняется, цель снова 100%
                // с < 0.99.. - немного не добрали в прошлый раз, но цель была добрать максимум, поэтому не меняем ничего - это погрешность
                // c > 2.1 - запасов хватит на три дня или более, нет никакой потребности
                y_int[i] = Math.max(0, n - y_int[i]);
            }
        }
        return y_int;
    }


}
