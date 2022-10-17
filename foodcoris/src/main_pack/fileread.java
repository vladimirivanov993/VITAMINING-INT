package main_pack;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.Vector;

public class fileread {

    public fileread() {
        // TODO Auto-generated constructor stub


    }

    public static void main(String[] args) throws FileNotFoundException {
        Scanner sc = new Scanner(new File("in.txt"));
        String[] ms;
        Vector<String> names = new Vector<String>();
        while (sc.hasNext()) {
            String s = sc.nextLine();
            ms = s.split(";");
            names.add(ms[0]);
            for (int i = 2; i < ms.length; i++) {

            }
        }
        sc.close();
        PrintWriter pw = new PrintWriter(new File("out.txt"));
        pw.close();


    }

}
