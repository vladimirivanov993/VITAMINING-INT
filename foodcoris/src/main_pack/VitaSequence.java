package main_pack;

import java.util.Vector;

public class VitaSequence {
    Vector<String> vs2 = new Vector<String>();

    double[] b; //стоимость продуктов
//	double[] b = {23,100,72,439,200,19,52,20,15,17,24,250,40,35,75,235,95,110,520,537,50,460,115,450,162,60,75,60,48,180,65,79,440,45,65,75,103,350,260,168,90,165,200,51,98,85,180,54,65,59,75,40,170,173,119,400,43,26};

    double[] c; //потребности и ограничения
//    double[] c = {100,80,65,320,-6000,75,0,1.5,1.5,2.0,20,2,0,0.4,0.25,0.2,0,20,1500,2000,1125,400,1000,19,0.1};

    Vector<String> vs; //линия строки данных - показатели количества витамин в продуктах

    Vector<String> names; //имена продуктов, попадающих в результат, количества витамин в результате; также это все строки, которые пойдут в основной отчет

    double[] massive_overfit; //вектор перегрузки значения

    public VitaSequence() {

        // TODO Auto-generated constructor stub


    }

}
